package zfactor

import "fmt"

// InputError represents an error resulting from invalid input parameters.
type InputError struct {
	Msg string
}

func (e InputError) Error() string {
	return e.Msg
}

// RangeError reports that a temperature lies outside the interval over
// which a correlation was fitted.
//
// The correlations in this module are fits to experimental data and stay
// mathematically defined beyond the interval they were fitted over, so
// the calculation is still performed and the extrapolated value returned
// alongside this error rather than discarded. Callers that accept the
// extrapolation may disregard it; those needing a value backed by data
// must check for it. The vapor-liquid equilibrium solvers take the
// former view, treating a RangeError as a caveat and any other error as
// fatal.
//
// Accuracy falls away with distance from the fitted interval, and a
// correlation is not constrained to behave sensibly far outside it.
//
// The type is shared so that a caller can test for an out-of-range
// result uniformly, whichever correlation produced it:
//
//	var rangeErr *zfactor.RangeError
//	if errors.As(err, &rangeErr) { ... }
type RangeError struct {
	Name string  // the substance or correlation the range belongs to
	T    float64 // the offending temperature
	Low  float64 // lower bound of the fitted range
	High float64 // upper bound of the fitted range
}

func (e *RangeError) Error() string {
	if e.Name == "" {
		return fmt.Sprintf(
			"temperature %g is outside the fitted range [%g, %g]",
			e.T, e.Low, e.High,
		)
	}

	return fmt.Sprintf(
		"temperature %g is outside the range [%g, %g] fitted for %s",
		e.T, e.Low, e.High, e.Name,
	)
}

var (
	// ErrTemp is returned when the absolute temperature is less than or equal to 0.
	ErrTemp = InputError{Msg: "absolute temperature (T) cannot be less than or equal to 0"}
	// ErrPressure is returned when the pressure is less than 0.
	ErrPressure = InputError{Msg: "pressure (P) cannot be less than 0"}
	// ErrCriticalProp is returned when a critical property (Tc or Pc) is less than or equal to 0.
	ErrCriticalProp = InputError{Msg: "critical property (Tc, Pc, Vc or Zc) cannot have a value less than or equal to 0"}
	// ErrUniversalConst is returned when the universal gas constant (R) is less than or equal to 0.
	ErrUniversalConst = InputError{Msg: "universal gas constant (R) value cannot be less than or equal to 0"}
	// ErrVirialCoeff is returned when a virial coefficient is 0.
	ErrVirialCoeff = InputError{Msg: "virial coefficient (B or C) cannot be 0"}
	// ErrVolume is returned when the molar volume is less than or equal to 0
	ErrVolume = InputError{Msg: "molar volume (V) cannot be less than or equal to 0"}
	// ErrHighPressureTwoTerm is returned when the pressure exceeds 15 bar for the two-term virial equation.
	ErrHighPressureTwoTerm = InputError{Msg: "pressure exceeds the validity limit (15 bar) for the two-term virial equation"}
	// ErrInvalidTr is returned when the reduced temperature (Tr) is less than or equal to 0.
	ErrInvalidTr = InputError{Msg: "reduced temperature (Tr) must be greater than 0"}
	// ErrInvalidPr is returned when the reduced pressure (Pr) is less than or equal to 0.
	ErrInvalidPr = InputError{Msg: "reduced pressure (Pr) must be greater than 0"}
	// ErrMolFracSum is returned when the mole fractions do not add up to 1 or are at least out of the tolerance range.
	ErrMolFracSum = InputError{Msg: "mole fractions should sum to 1.0"}
	// ErrMolFracVal is returned when the mole fraction is out of range.
	ErrMolFracVal = InputError{Msg: "mole fractions should sum to 1.0"}
)
