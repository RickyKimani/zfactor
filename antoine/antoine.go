// Package antoine provides coefficients and calculation methods for the Antoine equation,
// which estimates the saturation vapor pressure of pure substances as a function of temperature.
//
// The form used is: ln(P[kPa]) = A - B / (T[°C] + C)
package antoine

import (
	"fmt"
	"math"

	"github.com/rickykimani/zfactor"
)

// RangeError reports that a temperature lies outside the interval over
// which a correlation was fitted. It is an alias for zfactor.RangeError,
// so a caller may test for an out-of-range result from any correlation
// in this module with a single type.
//
// Temperatures here are in °C. The equation diverges at t = -C, which
// lies below the fitted range of every tabulated substance.
type RangeError = zfactor.RangeError

// Model evaluates the saturation vapor pressure of a pure substance as a
// function of temperature, together with the inverse relation.
//
// The interface exists so that vapor-liquid equilibrium calculations can
// accept any Antoine-like correlation without knowing which. Antoine is
// the implementation provided here; a caller may supply an extended
// Antoine form, a Wagner equation or a fit of their own in its place.
//
// Implementations work in degrees Celsius and kilopascals, matching the
// coefficients tabulated in this package.
//
// A temperature outside the fitted range is reported with a *RangeError
// returned alongside the computed value rather than in place of it, so
// that callers willing to extrapolate can proceed. Any other error means
// no value could be produced.
type Model interface {
	// LnPSat returns the natural logarithm of the saturation pressure
	// in kPa at temperature t, in °C.
	LnPSat(t float64) (float64, error)

	// Pressure returns the saturation pressure in kPa at temperature t,
	// in °C.
	Pressure(t float64) (float64, error)

	// ValidateTempRange reports whether t, in °C, lies within the
	// interval the correlation was fitted over. Both bounds are
	// inclusive.
	ValidateTempRange(t float64) bool

	// Temperature returns the saturation temperature in °C at pressure
	// p, in kPa, inverting the correlation.
	//
	// It returns an error when no saturation temperature corresponds to
	// p, which is the case for a non-positive pressure and for one large
	// enough that the inverted equation has no solution.
	Temperature(p float64) (float64, error)
}

// Antoine holds the constants for the Antoine equation: ln(P) = A - B/(T+C)
// Units: P in kPa, T in °C
type Antoine struct {
	Name    string
	Formula string
	A       float64
	B       float64
	C       float64
	H       float64   // Latent heat of vaporization (kJ/mol)
	Range   TempRange // Valid temperature range (°C)
	Tn      float64   // Normal boiling point (°C)
}

// TempRange defines a valid temperature interval.
type TempRange struct {
	Low  float64
	High float64
}

// LnPSat calculates the natural logarithm of the saturation pressure (kPa) at temperature t (°C).
// Returns an error if t is outside the valid range.
func (a *Antoine) LnPSat(t float64) (float64, error) {
	var err error
	if !a.ValidateTempRange(t) {
		err = &RangeError{
			Name: a.Name,
			T:    t,
			Low:  a.Range.Low,
			High: a.Range.High,
		}
	}
	return a.A - a.B/(t+a.C), err
}

// Pressure calculates the saturation pressure (kPa) at temperature t (°C).
// Returns an error if t is outside the valid range.
func (a *Antoine) Pressure(t float64) (float64, error) {
	lnP, err := a.LnPSat(t)

	return math.Exp(lnP), err
}

// ValidateTempRange reports whether t lies within the valid temperature range.
func (a *Antoine) ValidateTempRange(t float64) bool {
	return t >= a.Range.Low && t <= a.Range.High
}

// Temperature calculates the saturation temperature (°C) at a pressure p (kPa).
//
// Inverting ln(P) = A - B/(t + C) gives t = B/(A - ln P) - C, which has
// no solution once ln P reaches A: the denominator vanishes there and
// changes sign beyond it, so the expression returns an infinity or a
// temperature below the pole at -C. Both are rejected rather than
// returned, since neither describes a saturation state.
//
// Returns an error if p is non-positive or at or beyond that limit.
func (a *Antoine) Temperature(p float64) (float64, error) {
	if p <= 0 {
		return 0, zfactor.ErrPressure
	}

	denominator := a.A - math.Log(p)
	if denominator <= 0 {
		return 0, fmt.Errorf(
			"pressure %g kPa is at or beyond the limit of the correlation for %s (ln p must stay below A = %g)",
			p, a.Name, a.A,
		)
	}

	return a.B/denominator - a.C, nil
}
