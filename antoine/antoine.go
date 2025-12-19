package antoine

import (
	"fmt"
	"math"
)

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
	if t < a.Range.Low || t > a.Range.High {
		return 0, fmt.Errorf("temperature %.2f°C is outside the valid range [%.2f, %.2f]", t, a.Range.Low, a.Range.High)
	}
	return a.A - a.B/(t+a.C), nil
}

// Pressure calculates the saturation pressure (kPa) at temperature t (°C).
// Returns an error if t is outside the valid range.
func (a *Antoine) Pressure(t float64) (float64, error) {
	lnP, err := a.LnPSat(t)
	if err != nil {
		return 0, err
	}
	return math.Exp(lnP), nil
}
