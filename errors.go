package zfactor

import "errors"

var (
	// TempErr is returned when the absolute temperature is less than or equal to 0.
	TempErr = errors.New("absolute temperature (T) cannot be less than or equal to 0")
	// PressErr is returned when the pressure is less than 0.
	PressErr = errors.New("pressure (P) cannot be less than 0")
	// CriticalPropErr is returned when a critical property (Tc or Pc) is less than or equal to 0.
	CriticalPropErr = errors.New("critical property (Tc or Pc) cannot have a value less than or equal to 0")
	// UniversalConstErr is returned when the universal gas constant (R) is less than or equal to 0.
	UniversalConstErr = errors.New("universal gas constant (R) value cannot be less than or equal to 0")
	// VirialCoeffErr is returned when a virial coefficient is 0.
	VirialCoeffErr = errors.New("virial coefficient (B or C) cannot be 0")
	// VolumeErr is returned when the molar volume is less than or equal to 0
	VolumeErr      = errors.New("molar volume (V) cannot be less than or equal to 0")
)
