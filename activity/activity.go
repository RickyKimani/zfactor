// Package activity provides thermodynamic activity coefficient models and
// utilities for liquid-phase solution thermodynamics.
//
// The package is organized into model-specific subpackages implementing
// common excess Gibbs energy models, such as:
//
//   - Wilson
//   - NRTL
//   - Margules
//   - van Laar
//
// Each model exposes:
//
//   - A Data type containing the model parameters and system state.
//   - An Activity function for evaluating activity coefficients.
//
// Generic numerical utilities that are independent of the underlying
// thermodynamic model, such as infinite-dilution calculations, are provided
// by the root activity package.
//
// Example:
//
//	data := wilson.Data{
//	    X:           x,
//	    Interaction: lambda,
//	    V:           v,
//	    T:           T,
//	}
//
//	gamma := func(x []float64) ([]float64, error) {
//	    d := data
//	    d.X = x
//	    return wilson.Activity(d)
//	}
//
//	gammaInf, err := activity.InfiniteDilution(data.X, gamma)
package activity

import (
	"fmt"
)

const Tolerance = 1e-9

// Model represents a liquid-phase activity coefficient model.
//
// Implementations encapsulate the parameters required to evaluate
// activity coefficients for a particular excess Gibbs energy model
// (e.g. Wilson or NRTL). The interface is designed so that generic
// algorithms, such as infinite-dilution calculations, can operate on
// any supported model without knowledge of its underlying formulation.
type Model interface {
	// Activity returns the activity coefficients for the model's
	// current liquid-phase composition.
	Activity() ([]float64, error)

	// Composition returns the current liquid-phase mole fraction vector.
	//
	// Implementations should return a copy of the underlying slice to
	// prevent callers from modifying the model's state.
	Composition() []float64

	// WithComposition returns a copy of the model with the supplied
	// liquid-phase composition.
	//
	// The receiver must not be modified.
	WithComposition([]float64) Model
}

// BinaryInfiniteDilutioner is implemented by models that provide an
// analytical closed-form solution for binary infinite-dilution activity
// coefficients.
//
// When available, generic algorithms may use this method in preference
// to numerical approximation.
type BinaryInfiniteDilutioner interface {
	// BinaryInfiniteDilution returns the infinite-dilution activity
	// coefficients for a binary mixture.
	//
	// The returned slice contains γ∞ values in the same order as the
	// model's composition vector.
	BinaryInfiniteDilution() ([]float64, error)
}

// InfiniteDilution approximates the infinite-dilution activity coefficients
// for every component in a mixture.
//
// For each component i, the liquid composition is perturbed such that
//
//	x_i = ε
//
// while the remaining mole fractions are scaled proportionally so that
//
//	Σ x = 1.
//
// A new model with the perturbed composition is created using
// Model.WithComposition, and the activity coefficients are evaluated.
// The activity coefficient of the diluted component is taken as its
// infinite-dilution value.
//
// The supplied model is never modified.
func InfiniteDilution(model Model) ([]float64, error) {
	const eps = 1e-12

	x := model.Composition()
	m := len(x)

	if m == 0 {
		return nil, fmt.Errorf("no components provided")
	}

	if b, ok := model.(BinaryInfiniteDilutioner); ok && m == 2 {
		return b.BinaryInfiniteDilution()
	}

	gammaInf := make([]float64, m)

	for i := range m {
		xc := make([]float64, m)
		copy(xc, x)

		rem := 0.0
		for j := range m {
			if j != i {
				rem += xc[j]
			}
		}

		if rem == 0 {
			return nil, fmt.Errorf(
				"cannot dilute component %d: remaining composition is zero",
				i,
			)
		}

		// Dilute component i.
		xc[i] = eps

		// Renormalize remaining components.
		scale := (1 - eps) / rem

		for j := range m {
			if j != i {
				xc[j] *= scale
			}
		}

		gamma, err := model.WithComposition(xc).Activity()
		if err != nil {
			return nil, err
		}

		if len(gamma) != m {
			return nil, fmt.Errorf(
				"activity returned %d coefficients, expected %d",
				len(gamma),
				m,
			)
		}

		gammaInf[i] = gamma[i]
	}

	return gammaInf, nil
}
