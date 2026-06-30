// Package nrtl implements the Non-Random Two-Liquid (NRTL)
// activity coefficient model.
//
// The NRTL model, proposed by Renon and Prausnitz (1968),
// is a local-composition model for representing non-ideal
// liquid mixtures. It is widely used for vapor-liquid and
// liquid-liquid equilibrium calculations.
//
// The model calculates activity coefficients from binary
// interaction parameters (τ) and non-randomness parameters
// (α). Temperature-dependent interaction parameters may be
// evaluated using the extended correlation:
//
//	τᵢⱼ = aᵢⱼ + bᵢⱼ/T + cᵢⱼ ln(T) + dᵢⱼT
package nrtl

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/activity"
)

const tol = activity.Tolerance

// TauCorrelation represents a source of NRTL binary interaction
// parameters.
//
// Implementations may return constant interaction parameters or
// evaluate temperature-dependent correlations.
//
// The supplied temperature must be in Kelvin.
type TauCorrelation interface {
    Tau(T float64) ([][]float64, error)
}

// ConstantTau represents temperature-independent binary interaction
// parameters.
type ConstantTau struct {
    TauMatrix [][]float64
}

// Tau returns the stored interaction parameter matrix.
// The temperature argument is ignored.
func (c ConstantTau) Tau(float64) ([][]float64, error) {
    return c.TauMatrix, nil
}

// ExtendedTau evaluates NRTL interaction parameters
// from the extended temperature-dependent correlation
//
//    τij = aij + bij/T + cij ln(T) + dijT
//
// where T is the absolute temperature in Kelvin.
type ExtendedTau struct {
	A [][]float64
	B [][]float64
	C [][]float64
	D [][]float64
}

// Tau evaluates the NRTL interaction parameter matrix at the
// specified absolute temperature T (K).
func (td ExtendedTau) Tau(T float64) ([][]float64, error) {
	if T <= 0 {
		return nil, zfactor.ErrTemp
	}

	a := td.A
	b := td.B
	c := td.C
	d := td.D

	m := len(a)
	if m == 0 {
		return nil, errors.New("no components provided")
	}
	if len(b) != m {
		return nil, errors.New("incorrect number of b parameter rows provided")
	}
	if len(c) != m {
		return nil, errors.New("incorrect number of c parameter rows provided")
	}
	if len(d) != m {
		return nil, errors.New("incorrect number of d parameter rows provided")
	}
	for i := range m {
		if len(a[i]) != m {
			return nil, errors.New("incorrect number of a parameters provided")
		}
		if len(b[i]) != m {
			return nil, errors.New("incorrect number of b parameters provided")
		}
		if len(c[i]) != m {
			return nil, errors.New("incorrect number of c parameters provided")
		}
		if len(d[i]) != m {
			return nil, errors.New("incorrect number of d parameters provided")
		}
	}

	tau := make([][]float64, m)
	for i := range m {
		tau[i] = make([]float64, m)
		for j := range m {
			tau[i][j] = a[i][j] + b[i][j]/T + c[i][j]*math.Log(T) + d[i][j]*T
		}
	}

	return tau, nil
}

// NRTL represents an N-component liquid mixture described by the
// Non-Random Two-Liquid model.
//
// Temperature is specified in degrees Celsius. Temperature-dependent
// interaction parameters are obtained through the configured
// TauCorrelation implementation.
type NRTL struct {

	// System temperature
	T float64

	// Mole fractions of each component.
	X []float64

	// Non-randomness parameter matrix (αij).
	Alpha [][]float64

	// Binary interaction parameter matrix (τij).
	Tau TauCorrelation
}

// Composition returns a copy of the liquid-phase mole fraction vector.
func (n NRTL) Composition() []float64 {
	x := make([]float64, len(n.X))
	copy(x, n.X)
	return x
}

// Temperature returns the temperature supplied to the model
func (n NRTL) Temperature() float64 {
	return n.T
}

// WithComposition returns a copy of the NRTL model with the supplied
// liquid-phase composition.
func (n NRTL) WithComposition(x []float64) activity.Model {
	n.X = make([]float64, len(x))
	copy(n.X, x)
	return n
}

// WithTemperature returns a copy of the model with the supplied temperature.
func (n NRTL) WithTemperature(T float64) activity.Model {
	n.T = T
	return n
}

// tau evaluates the interaction parameter matrix at the model
// temperature.
func (n NRTL) tau() ([][]float64, error) {
    if n.Tau == nil {
        return nil, errors.New("no tau model provided")
    }

    tau, err := n.Tau.Tau(n.T)
    if err != nil {
        return nil, err
    }

    return tau, nil
}

// Activity calculates the activity coefficients of all
// components in a liquid mixture using the NRTL model.
//
// The mole fractions must sum to unity within the package
// tolerance. Alpha and Tau must both be square n×n matrices,
// where m is the number of components.
//
// The returned slice contains the activity coefficient of
// each component in the same order as the input composition.
func (n NRTL) Activity() ([]float64, error) {
	m := len(n.X)
	if m == 0 {
		return nil, errors.New("no components provided")
	}

	if len(n.Alpha) != m {
		return nil, errors.New("incorrect number of alpha rows")
	}

	tau, err := n.tau()
	if err != nil {
		return nil, err
	}

	if len(tau) != m {
		return nil, errors.New("incorrect number of tau rows")
	}

	for i := range m {
		if len(n.Alpha[i]) != m {
			return nil, errors.New("incorrect number of alpha parameters")
		}
		if len(tau[i]) != m {
			return nil, errors.New("incorrect number of tau parameters")
		}
	}

	x := n.X
	alpha := n.Alpha

	// Check mole fractions
	sumX := 0.0
	for _, val := range x {
		if val < 0 || val > 1 {
			return nil, zfactor.ErrMolFracVal
		}
		sumX += val
	}
	if math.Abs(sumX-1.0) > tol {
		return nil, zfactor.ErrMolFracSum
	}

	// Compute the NRTL weighting factors:
	//
	//	Gij = exp(-αijτij)
	G := make([][]float64, m)
	for i := range m {
		G[i] = make([]float64, m)
		for j := range m {
			G[i][j] = math.Exp(-alpha[i][j] * tau[i][j])
		}
	}

	// Compute
	//
	//	Σk xk Gki
	//
	// for each component.
	weightedG := make([]float64, m)
	for i := range m {
		for k := range m {
			weightedG[i] += x[k] * G[k][i]
		}

		if weightedG[i] == 0 {
			return nil, errors.New("division by zero in GSum")
		}
	}

	// Compute
	//
	//	Σj xjτjiGji
	//
	// for each component.
	weightedTauG := make([]float64, m)
	for i := range m {
		for j := range m {
			weightedTauG[i] += x[j] * tau[j][i] * G[j][i]
		}
	}

	// Evaluate the NRTL activity coefficients.
	gamma := make([]float64, m)

	for i := range m {

		// First contribution:
		//
		//	Σj xjτjiGji
		//	────────────
		//	  Σk xkGki
		term1 := weightedTauG[i] / weightedG[i]

		// Second contribution:
		//
		//	    xjGij
		//	────────────── x (τij - Σm xmτmjGmj / Σk xkGkj)
		//	  Σk xkGkj
		term2 := 0.0

		for j := range m {

			prefactor := x[j] * G[i][j] / weightedG[j]

			bracket := tau[i][j] - weightedTauG[j]/weightedG[j]

			term2 += prefactor * bracket
		}

		lnGamma := term1 + term2

		gamma[i] = math.Exp(lnGamma)
	}

	return gamma, nil
}

// BinaryInfiniteDilution returns the closed-form infinite-dilution
// activity coefficients for a binary mixture using the NRTL model.
//
// For a binary system,
//
//	G12 = exp(-α12τ12)
//	G21 = exp(-α21τ21)
//
// and
//
//	ln(γ1∞) = τ21 + G12(τ12 - τ21)
//	ln(γ2∞) = τ12 + G21(τ21 - τ12)
//
// The receiver must contain exactly two components.
func (n NRTL) BinaryInfiniteDilution() ([]float64, error) {
	m := 2
	if len(n.X) != m {
		return nil, errors.New("binary infinite dilution requires exactly two components")
	}

	if len(n.Alpha) != m {
		return nil, errors.New("incorrect number of alpha rows")
	}

	tau, err := n.tau()
	if err != nil {
		return nil, err
	}

	if len(tau) != m {
		return nil, errors.New("incorrect number of tau rows")
	}

	for i := range m {
		if len(n.Alpha[i]) != m {
			return nil, errors.New("incorrect number of alpha parameters")
		}
		if len(tau[i]) != m {
			return nil, errors.New("incorrect number of tau parameters")
		}
	}

	g12 := math.Exp(-n.Alpha[0][1] * tau[0][1])
	g21 := math.Exp(-n.Alpha[1][0] * tau[1][0])

	lnGamma1 := tau[1][0] +
		g12*(tau[0][1]-tau[1][0])

	lnGamma2 := tau[0][1] +
		g21*(tau[1][0]-tau[0][1])

	return []float64{
		math.Exp(lnGamma1),
		math.Exp(lnGamma2),
	}, nil
}
