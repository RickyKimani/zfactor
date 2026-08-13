// Package cp provides coefficients for the calculation of heat capacities
// of gases (in the ideal-gas state), solids, and liquids, together with
// the ideal-gas enthalpy and entropy changes obtained by integrating them.
//
// The parameters stored in HeatCapacity struct correspond to the equation:
//
//	Cp/R = A + B*T + C*T^2 + D*T^-2
package cp

import (
	"errors"
	"fmt"
	"math"

	"github.com/rickykimani/zfactor"
)

// HeatCapacity holds the constants for the heat capacity equation.
// Cp/R = A + B*T + C*T^2 + D*T^-2
type HeatCapacity struct {
	Name    string
	Formula string
	TMin    float64
	TMax    float64
	Cp298   float64 // Cp/R at 298.15 K
	A       float64
	B       float64
	C       float64
	D       float64
}

// RangeError reports that a temperature lies outside the interval over
// which a heat-capacity correlation was fitted. It is an alias for
// zfactor.RangeError, so a caller may test for an out-of-range result
// from any correlation in this module with a single type.
//
// Temperatures here are in kelvin.
type RangeError = zfactor.RangeError

// checkRange reports whether T lies within the fitted interval.
func (h *HeatCapacity) checkRange(T float64) error {
	if T < h.TMin || T > h.TMax {
		return &RangeError{Name: h.Name, T: T, Low: h.TMin, High: h.TMax}
	}
	return nil
}

// IdealGasEnthalpyChange calculates the change in enthalpy (Delta H) for an ideal gas state
// between two states.
//
// The state arguments must include:
//   - T: Temperature (K)
//   - P: Pressure
//   - R: Universal Gas Constant
//
// Formula: Delta H = R * Integral(Cp/R dT) from T1 to T2
//
// A temperature outside the correlation's fitted range does not prevent
// the calculation: the extrapolated value is returned together with a
// *RangeError, which callers may inspect with errors.As. Invalid input —
// a non-positive temperature or pressure, or conflicting gas constants —
// yields zero and an error instead, since no meaningful value exists.
func (h *HeatCapacity) IdealGasEnthalpyChange(state1, state2 zfactor.Args) (float64, error) {
	T1 := state1.T
	P1 := state1.P

	T2 := state2.T
	P2 := state2.P

	if state1.R != state2.R {
		return 0, fmt.Errorf("conflicting values of R {%v and %v}", state1.R, state2.R)
	}

	R := state1.R

	if T1 <= 0 || T2 <= 0 {
		return 0, zfactor.ErrTemp
	}
	if P1 <= 0 || P2 <= 0 {
		return 0, zfactor.ErrPressure
	}

	// Reported alongside the result rather than in place of it.
	rangeErr := errors.Join(h.checkRange(T1), h.checkRange(T2))

	// Integral of Cp/R * dT
	// Cp/R = A + BT + CT^2 + DT^-2
	// Int  = AT + (B/2)T^2 + (C/3)T^3 - D/T

	termA := h.A * (T2 - T1)
	termB := (h.B / 2) * (T2*T2 - T1*T1)
	termC := (h.C / 3) * (T2*T2*T2 - T1*T1*T1)
	termD := -h.D * (1/T2 - 1/T1)

	return R * (termA + termB + termC + termD), rangeErr
}

// IdealGasEntropyChange calculates the change in entropy (Delta S) for an ideal gas state
// between two states.
//
// The state arguments must include:
//   - T: Temperature (K)
//   - P: Pressure
//   - R: Universal Gas Constant
//
// Formula: Delta S = R * (Integral(Cp/(R*T) dT) - ln(P2/P1))
//
// As with IdealGasEnthalpyChange, a temperature outside the fitted range
// returns the extrapolated value together with a *RangeError, while
// invalid input returns zero and an error.
func (h *HeatCapacity) IdealGasEntropyChange(state1, state2 zfactor.Args) (float64, error) {
	T1 := state1.T
	P1 := state1.P

	T2 := state2.T
	P2 := state2.P

	if state1.R != state2.R {
		return 0, fmt.Errorf("conflicting values of R {%v and %v}", state1.R, state2.R)
	}

	R := state1.R

	if T1 <= 0 || T2 <= 0 {
		return 0, zfactor.ErrTemp
	}
	if P1 <= 0 || P2 <= 0 {
		return 0, zfactor.ErrPressure
	}

	// Reported alongside the result rather than in place of it.
	rangeErr := errors.Join(h.checkRange(T1), h.checkRange(T2))

	// Integral of (Cp/R)/T * dT
	// Cp/R = A + BT + CT^2 + DT^-2
	// (Cp/R)/T = A/T + B + CT + DT^-3
	// Int = A*ln(T) + BT + (C/2)T^2 - (D/2)T^-2

	termA := h.A * (math.Log(T2 / T1))
	termB := h.B * (T2 - T1)
	termC := (h.C / 2) * (T2*T2 - T1*T1)
	// For D term: - (D/2) * (1/T2^2 - 1/T1^2)
	termD := -(h.D / 2) * ((1 / (T2 * T2)) - (1 / (T1 * T1)))

	integral := termA + termB + termC + termD
	pressureTerm := math.Log(P2 / P1)

	return R * (integral - pressureTerm), rangeErr
}
