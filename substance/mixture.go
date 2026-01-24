package substance

import (
	"errors"
	"math"
)

// Component represents a pure substance and its mole fraction in a mixture.
type Component struct {
	Substance *Substance
	Fraction  float64
}

// NewLinearMixture creates a new 'pseudo-substance' representing a mixture of gases.
// This is mainly intended for gas mixtures to estimate properties using generalized correlations.
//
// It uses Kay's Rule (linear molar averages) to calculate pseudo-critical properties,
// molar mass, and acentric factor.
//
//	Tc_mix = Σ (yi * Tci)
//	Pc_mix = Σ (yi * Pci)
//	ω_mix  = Σ (yi * ωi)
//	MW_mix = Σ (yi * MWi)
//	Vc_mix = Σ (yi * Vci)
//	Zc_mix = Σ (yi * Zci)
//
// This approach treats the mixture as a single fictional "pseudo-component". This is
// widely used with generalized correlations like Lee-Kesler.
func NewLinearMixture(name string, components []Component) (*Substance, error) {
	if len(components) == 0 {
		return nil, errors.New("mixture must have at least one component")
	}

	var (
		sumF     float64
		mix      Substance
		critical CriticalProps
	)
	mix.Name = name
	// Normal boiling point is not applicable to linear mixtures, so we set it to -1
	// to indicate it is unusable.
	mix.Tn = -1

	for _, c := range components {
		if c.Substance == nil {
			return nil, errors.New("component substance cannot be nil")
		}
		if c.Fraction < 0 {
			return nil, errors.New("mole fraction cannot be negative")
		}

		y := c.Fraction
		sumF += y

		// Linear averages
		mix.MW += y * c.Substance.MW
		mix.Acentric += y * c.Substance.Acentric

		critical.Tc += y * c.Substance.Critical.Tc
		critical.Pc += y * c.Substance.Critical.Pc
		critical.Vc += y * c.Substance.Critical.Vc
		critical.Zc += y * c.Substance.Critical.Zc
	}

	const tolerance = 1e-4
	if math.Abs(sumF-1.0) > tolerance {
		return nil, errors.New("mole fractions must sum to 1.0")
	}

	mix.Critical = critical
	return &mix, nil
}
