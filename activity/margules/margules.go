// Package margules implements the two-parameter Margules activity
// coefficient model for binary liquid mixtures.
package margules

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/activity"
)

// Margules represents the two-parameter Margules model for a binary
// liquid mixture.
//
// A12 and A21 are the Margules interaction parameters.
//
// X contains the liquid-phase mole fractions [x1, x2].
type Margules struct {
	A12 float64
	A21 float64
	X   []float64
}

// Activity calculates the liquid-phase activity coefficients using the
// two-parameter Margules model.
//
// For a binary mixture,
//
//	ln(γ₁) = x₂²[A₁₂ + 2(A₂₁−A₁₂)x₁]
//
//	ln(γ₂) = x₁²[A₂₁ + 2(A₁₂−A₂₁)x₂]
//
// The returned slice contains the activity coefficients
// [γ₁, γ₂].
func (m Margules) Activity() ([]float64, error) {
	if len(m.X) != 2 {
		return nil, errors.New("margules model requires exactly two components")
	}

	sum := 0.0
	for _, xi := range m.X {
		if xi < 0 || xi > 1 {
			return nil, zfactor.ErrMolFracVal
		}
		sum += xi
	}

	if math.Abs(sum-1) > activity.Tolerance {
		return nil, zfactor.ErrMolFracSum
	}

	x1 := m.X[0]
	x2 := m.X[1]

	lnGamma1 := x2 * x2 * (m.A12 + 2*(m.A21-m.A12)*x1)
	lnGamma2 := x1 * x1 * (m.A21 + 2*(m.A12-m.A21)*x2)

	return []float64{
		math.Exp(lnGamma1),
		math.Exp(lnGamma2),
	}, nil
}

// Composition returns a copy of the liquid-phase mole fraction vector.
func (m Margules) Composition() []float64 {
	x := make([]float64, len(m.X))
	copy(x, m.X)
	return x
}

// WithComposition returns a copy of the model with the supplied
// liquid-phase composition.
func (m Margules) WithComposition(x []float64) activity.Model {
	m.X = make([]float64, len(x))
	copy(m.X, x)
	return m
}

// BinaryInfiniteDilution returns the infinite-dilution activity
// coefficients for the Margules model.
//
// The closed-form expressions are
//
//	ln(γ₁∞) = A₁₂
//	ln(γ₂∞) = A₂₁
//
// For the symmetric one-parameter Margules model,
//
//	A₁₂ = A₂₁ = A₀
//
// giving
//
//	γ₁∞ = γ₂∞ = exp(A₀).
func (m Margules) BinaryInfiniteDilution() ([]float64, error) {
	return []float64{
		math.Exp(m.A12),
		math.Exp(m.A21),
	}, nil
}
