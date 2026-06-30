// Package wilson implements the Wilson activity coefficient model for
// predicting liquid-phase non-ideal behavior.
//
// The Wilson model expresses the excess Gibbs free energy using binary
// interaction parameters and pure-component liquid molar volumes. It is
// suitable for completely miscible liquid mixtures and provides activity
// coefficients for vapor-liquid equilibrium and related calculations.
//
// The primary entry points are Activity, which evaluates activity
// coefficients for a given mixture composition, and WilsonInfiniteDilution,
// which computes closed-form infinite-dilution activity coefficients.
package wilson

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/activity"
)

const tol = activity.Tolerance

// Wilson contains the data required to calculate Wilson liquid
// activity coefficients for every component in a multicomponent mixture.
//
// T is the system temperature.
//
// X is the liquid-phase mole fraction vector for the mixture.
//
// V is the liquid molar volume vector for the mixture.
//
// Interaction is the Wilson interaction parameter matrix aij (J/mol) such
// that Interaction[i][j] gives the parameter for component i with respect
// to component j.
//
// For an n-component mixture:
//
//	Interaction = [][]float64{
//	    {a11, a12, ..., a1n},
//	    {a21, a22, ..., a2n},
//	    ...
//	    {an1, an2, ..., ann},
//	}
//
// The diagonal terms aii are ignored because:
//
//	Λii = 1
type Wilson struct {
	T           float64     // system temperature (K)
	X           []float64   // composition vector
	V           []float64   // molar volume vector
	Interaction [][]float64 // interaction matrix
}

// Composition returns a copy of the liquid-phase mole fraction vector.
func (w Wilson) Composition() []float64 {
	x := make([]float64, len(w.X))
	copy(x, w.X)
	return x
}

// Temperature returns the temperature supplied to the model
func (w Wilson) Temperature() float64 {
	return w.T
}

// WithComposition returns a copy of the Wilson model with the supplied
// liquid-phase composition.
func (w Wilson) WithComposition(x []float64) activity.Model {
	w.X = make([]float64, len(x))
	copy(w.X, x)
	return w
}

// WithTemperature returns a copy of the model with the supplied temperature.
func (w Wilson) WithTemperature(T float64) activity.Model {
	w.T = T
	return w
}

// Activity calculates liquid-phase activity coefficients using
// the Wilson excess Gibbs energy model.
//
// The Wilson model expresses the excess Gibbs free energy as:
//
//	Gᴱ/(RT) = -Σ xi ln(Σ xj Λij)
//
// where
//
//	Λij = 1,                                    for i = j
//	Λij = (Vj/Vi) exp(-aij/(RT)),               for i != j
//
// Activity coefficients are obtained from:
//
//	ln(γi) = 1
//	         - ln(Σ xj Λij)
//	         - Σ [xk Λki / (Σ xj Λkj)]
//
// The Wilson model is commonly used for strongly non-ideal liquid
// mixtures and assumes complete miscibility of all components.
//
// All components must:
//
//   - Have positive molar volumes.
//   - Have mole fractions summing to unity.
//   - Provide an n x n Wilson parameter matrix for an n-component mixture.
//
// The returned slice contains activity coefficients in the same order
// as the input composition vector.
func (w Wilson) Activity() ([]float64, error) {
	R := zfactor.RSI
	m := len(w.X)
	if m == 0 {
		return nil, errors.New("no components provided")
	}

	// Validate inputs
	if w.T <= 0 {
		return nil, zfactor.ErrTemp
	}
	if len(w.V) != m {
		return nil, errors.New("incorrect number of molar volumes")
	}
	if len(w.Interaction) != m {
		return nil, errors.New("incorrect number of wilson parameter rows")
	}
	for i := range m {
		if w.V[i] <= 0 {
			return nil, zfactor.ErrVolume
		}
		if len(w.Interaction[i]) != m {
			return nil, errors.New("incorrect number of wilson parameters")
		}
	}

	x := w.X
	v := w.V
	a := w.Interaction
	T := w.T

	// validate mole frac
	sumF := 0.0
	for _, val := range x {
		if val < 0 || val > 1 {
			return nil, zfactor.ErrMolFracVal
		}
		sumF += val
	}
	if math.Abs(sumF-1) > tol {
		return nil, zfactor.ErrMolFracSum
	}

	// Calculate Lambda matrix
	lambda := make([][]float64, m)
	for i := range m {
		lambda[i] = make([]float64, m)
		for j := range m {
			if i == j {
				lambda[i][j] = 1.0
			} else {
				lambda[i][j] = (v[j] / v[i]) * math.Exp(-a[i][j]/(R*T))
			}
		}
	}

	denom := make([]float64, m)
	for k := range m {
		for j := range m {
			denom[k] += x[j] * lambda[k][j]
		}
	}

	// Calculate activity coefficients
	gamma := make([]float64, m)
	for i := range m {
		sum1 := denom[i]

		sum2 := 0.0
		for k := range m {
			sum2 += (x[k] * lambda[k][i]) / denom[k]
		}

		lnGamma := 1 - math.Log(sum1) - sum2
		gamma[i] = math.Exp(lnGamma)
	}

	return gamma, nil
}

// BinaryInfiniteDilution returns the closed-form infinite-dilution
// activity coefficients for a binary mixture using the Wilson model.
//
// The Wilson model gives
//
//	ln(γ₁∞) = -ln(Λ₁₂) + 1 - Λ₂₁
//	ln(γ₂∞) = -ln(Λ₂₁) + 1 - Λ₁₂
//
// where
//
//	Λᵢⱼ = (Vⱼ / Vᵢ) exp(-aᵢⱼ / RT)
//
// The receiver must contain exactly two components.
func (w Wilson) BinaryInfiniteDilution() ([]float64, error) {
	m := 2
	if len(w.X) != m {
		return nil, errors.New("binary infinite dilution requires exactly two components")
	}

	// Validate inputs
	if w.T <= 0 {
		return nil, zfactor.ErrTemp
	}
	if len(w.V) != m {
		return nil, errors.New("incorrect number of molar volumes")
	}
	if len(w.Interaction) != m {
		return nil, errors.New("incorrect number of wilson parameter rows")
	}
	for i := range m {
		if w.V[i] <= 0 {
			return nil, zfactor.ErrVolume
		}
		if len(w.Interaction[i]) != m {
			return nil, errors.New("incorrect number of wilson parameters")
		}
	}
	R := zfactor.RSI

	lambda12 := (w.V[1] / w.V[0]) *
		math.Exp(-w.Interaction[0][1]/(R*w.T))

	lambda21 := (w.V[0] / w.V[1]) *
		math.Exp(-w.Interaction[1][0]/(R*w.T))

	lnGamma1 := -math.Log(lambda12) + 1 - lambda21
	lnGamma2 := -math.Log(lambda21) + 1 - lambda12

	return []float64{
		math.Exp(lnGamma1),
		math.Exp(lnGamma2),
	}, nil
}
