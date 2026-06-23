// Package vanlaar implements the two-parameter Van Laar activity
// coefficient model for binary liquid mixtures.
package vanlaar

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/activity"
)

// VanLaar represents the two-parameter Van Laar model for a binary
// liquid mixture.
//
// A12 and A21 are the Van Laar interaction parameters.
//
// X contains the liquid-phase mole fractions [x1, x2].
type VanLaar struct {
	A12 float64
	A21 float64
	X   []float64
}

// Activity calculates the liquid-phase activity coefficients using the
// two-parameter Van Laar model.
//
// For a binary mixture,
//
//	ln(γ₁) = A₁₂(1 + A₁₂x₁/(A₂₁x₂))⁻²
//
//	ln(γ₂) = A₂₁(1 + A₂₁x₂/(A₁₂x₁))⁻²
//
// The returned slice contains the activity coefficients
// [γ₁, γ₂].
func (v VanLaar) Activity() ([]float64, error) {
	if len(v.X) != 2 {
		return nil, errors.New("van laar model requires exactly two components")
	}

	sum := 0.0
	for _, xi := range v.X {
		if xi < 0 || xi > 1 {
			return nil, zfactor.ErrMolFracVal
		}
		sum += xi
	}

	if math.Abs(sum-1) > activity.Tolerance {
		return nil, zfactor.ErrMolFracSum
	}

	x1 := v.X[0]
	x2 := v.X[1]

	if x1 == 0 || x2 == 0 {
		return nil, errors.New("activity coefficients are undefined for pure-component compositions; use BinaryInfiniteDilution for infinite dilution")
	}

	if v.A12 == 0 || v.A21 == 0 {
		return nil, errors.New("van laar parameters must be non-zero")
	}

	lnGamma1 := v.A12 * math.Pow(
		1+(v.A12*x1)/(v.A21*x2),
		-2,
	)

	lnGamma2 := v.A21 * math.Pow(
		1+(v.A21*x2)/(v.A12*x1),
		-2,
	)

	return []float64{
		math.Exp(lnGamma1),
		math.Exp(lnGamma2),
	}, nil
}

// Composition returns a copy of the liquid-phase mole fraction vector.
func (v VanLaar) Composition() []float64 {
	x := make([]float64, len(v.X))
	copy(x, v.X)
	return x
}

// WithComposition returns a copy of the model with the supplied
// liquid-phase composition.
func (v VanLaar) WithComposition(x []float64) activity.Model {
	v.X = make([]float64, len(x))
	copy(v.X, x)
	return v
}

// BinaryInfiniteDilution returns the infinite-dilution activity
// coefficients for the Van Laar model.
//
// The closed-form expressions are
//
//	ln(γ₁∞) = A₁₂
//	ln(γ₂∞) = A₂₁
//
// For the symmetric one-parameter Van Laar model,
//
//	A₁₂ = A₂₁ = A₀
//
// giving
//
//	γ₁∞ = γ₂∞ = exp(A₀).
func (v VanLaar) BinaryInfiniteDilution() ([]float64, error) {
	return []float64{
		math.Exp(v.A12),
		math.Exp(v.A21),
	}, nil
}
