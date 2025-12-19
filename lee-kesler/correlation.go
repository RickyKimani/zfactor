package leekesler

// Property is a Lee-Kesler correlation family (Z, H, S, PHI).
type Property int

const (
	Z   Property = iota // Compressibility factor
	H                   // Residual enthalpy
	S                   // Residual entropy
	PHI                 // Fugacity coefficient
)

// correlation bundles the base ("0") and departure ("1") tables
// for a given property and exposes an At method to evaluate both.
type correlation struct {
	base   table // e.g., Z0, H0, S0, PHI0
	depart table // e.g., Z1, H1, S1, PHI1
}

// Correlation returns an evaluator for a property.
//
// Usage:
//
//	z0, z1, err := leekesler.Correlation(leekesler.Z).At(pr, tr)
func Correlation(p Property) correlation {
	switch p {
	case Z:
		return correlation{base: Z0Table, depart: Z1Table}
	case H:
		return correlation{base: H0Table, depart: H1Table}
	case S:
		return correlation{base: S0Table, depart: S1Table}
	case PHI:
		return correlation{base: PHI0Table, depart: PHI1Table}
	default:
		// Fallback to Z
		return correlation{base: Z0Table, depart: Z1Table}
	}
}

// At returns the base and departure values at (pr, tr).
// For Z, this returns (Z0, Z1).
func (c correlation) At(pr, tr float64) (float64, float64, error) {
	v0, err := c.base.At(pr, tr)
	if err != nil {
		return 0, 0, err
	}
	v1, err := c.depart.At(pr, tr)
	if err != nil {
		return 0, 0, err
	}
	return v0, v1, nil
}
