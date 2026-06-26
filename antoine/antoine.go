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

type RangeError struct {
	T    float64
	Low  float64
	High float64
}

func (r RangeError) Error() string {
	return fmt.Sprintf("t = %.2f is outside the range[%.2f-%.2f]", r.T, r.Low, r.High)
}

// TODO: Better doc for the interface

// Model is an interface all Antoine-like correlations implement
type Model interface {
	LnPSat(t float64) (float64, error)
	Pressure(t float64) (float64, error)
	ValidateTempRange(t float64) bool
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
// Returns an error if p is irregular.
func (a *Antoine) Temperature(p float64) (float64, error) {
	if p <= 0 {
		return 0, zfactor.ErrPressure
	}

	return a.B/(a.A-math.Log(p)) - a.C, nil
}
