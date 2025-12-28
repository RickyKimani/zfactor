package zfactor

// InputError represents an error resulting from invalid input parameters.
type InputError struct {
	Msg string
}

func (e InputError) Error() string {
	return e.Msg
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
)
