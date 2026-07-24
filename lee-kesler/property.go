package leekesler

import (
	"errors"
	"math"
)

// Property is a Lee-Kesler correlation family (Z, H, S, PHI). Each property
// bundles its base ("0") and departure ("1") tables together with the rule for
// folding those two values with the acentric factor (ω), so a caller never has
// to know which combination a given property uses.
type Property struct {
	base, depart *table
	combine      func(m0, m1, omega float64) float64
}

// additive folds the base and departure terms linearly: m0 + ω·m1.
// Used by Z, residual enthalpy, and residual entropy.
func additive(m0, m1, omega float64) float64 { return m0 + omega*m1 }

// powerLaw folds the base and departure terms as m0·m1^ω.
// Used by the fugacity coefficient.
func powerLaw(m0, m1, omega float64) float64 { return m0 * math.Pow(m1, omega) }

// The supported Lee-Kesler correlations.
var (
	CompressibilityFactor = Property{base: Z0Table, depart: Z1Table, combine: additive}     // Compressibility factor (Z)
	ResidualEnthalpy      = Property{base: H0Table, depart: H1Table, combine: additive}     // Dimensionless residual enthalpy (H^R / R*Tc)
	ResidualEntropy       = Property{base: S0Table, depart: S1Table, combine: additive}     // Dimensionless residual entropy (S^R / R)
	FugacityCoefficient   = Property{base: PHI0Table, depart: PHI1Table, combine: powerLaw} // Fugacity coefficient
)

// Eval evaluates the property at the given reduced temperature (Tr) and reduced
// pressure (Pr) for a fluid with the given acentric factor (omega). It looks up
// the base and departure table values and folds them with the property's rule.
//
// Usage:
//
//	z, err := leekesler.CompressibilityFactor.Eval(Tr, Pr, omega)
func (p Property) Eval(Tr, Pr, omega float64) (float64, error) {
	if p.base == nil || p.depart == nil || p.combine == nil {
		return 0, errors.New("leekesler: uninitialized Property")
	}

	m0, err := p.base.At(Tr, Pr)
	if err != nil {
		return 0, err
	}

	m1, err := p.depart.At(Tr, Pr)
	if err != nil {
		return 0, err
	}

	return p.combine(m0, m1, omega), nil
}
