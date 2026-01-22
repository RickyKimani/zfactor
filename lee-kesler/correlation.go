package leekesler

// Property is a Lee-Kesler correlation family (Z, H, S, PHI).
type Property int

const (
	CompressibilityFactor Property = iota // Compressibility factor (Z)
	ResidualEnthalpy                      // Dimensionless Residual enthalpy (H^R / R*Tc)
	ResidualEntropy                       // Dimensionless Residual entropy (S^R / R)
	FugacityCoefficient                   // Fugacity coefficient
)

// correlation bundles the base ("0") and departure ("1") tables
// for a given property and exposes an At method to evaluate both.
type correlation struct {
	base   *table // e.g., Z0, H0, S0, PHI0
	depart *table // e.g., Z1, H1, S1, PHI1
}

// Correlation returns an evaluator for a property.
//
// Usage:
//
//	z0, z1, err := leekesler.Correlation(leekesler.Z).At(Pr, Tr)
func Correlation(p Property) correlation {
	switch p {
	case CompressibilityFactor:
		return correlation{base: Z0Table, depart: Z1Table}
	case ResidualEnthalpy:
		return correlation{base: H0Table, depart: H1Table}
	case ResidualEntropy:
		return correlation{base: S0Table, depart: S1Table}
	case FugacityCoefficient:
		return correlation{base: PHI0Table, depart: PHI1Table}
	default:
		// Fallback to Z
		return correlation{base: Z0Table, depart: Z1Table} //panic instead?
	}
}

// At returns the base and departure values at (Tr, Pr).
// For Z, this returns (Z0, Z1).
func (c correlation) At(Tr, Pr float64) (float64, float64, error) {
	v0, err := c.base.At(Tr, Pr)
	if err != nil {
		return 0, 0, err
	}
	v1, err := c.depart.At(Tr, Pr)
	if err != nil {
		return 0, 0, err
	}
	return v0, v1, nil
}
