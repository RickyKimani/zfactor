// package substance contains the characteristic properties of pure
// species
package substance

import (
	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/abbott"
	"github.com/rickykimani/zfactor/cubic"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
	"github.com/rickykimani/zfactor/liquids"
)

type CriticalProps struct {
	Tc float64 //Critical Temperature (K)
	Pc float64 //Critical Pressure (bar)
	Vc float64 //Critical Volume (cm^3/mol)
	Zc float64 //Critical Compressibility factor
}

type Substance struct {
	Name     string
	MW       float64 //Molar mass
	Acentric float64 //Acentric factor
	Tn       float64 //Normal boiling point (K)
	Critical CriticalProps
}

// LeeKesler evaluates a thermodynamic property using the Lee-Kesler correlation.
// P is the pressure in bar.
// T is the temperature in Kelvin.
func (s *Substance) LeeKesler(T, P float64, property leekesler.Property) (float64, error) {
	pr := P / s.Critical.Pc
	tr := T / s.Critical.Tc

	c := leekesler.Correlation(property)

	v0, v1, err := c.At(tr, pr)
	if err != nil {
		return 0, err
	}

	v := v0 + s.Acentric*v1

	return v, nil
}

// CubicConfig creates a configuration for a cubic equation of state (EOS) solver.
// It initializes the EOS parameters based on the substance's critical properties and acentric factor.
//
// Supported standard types (VdW, RK, SRK, PR) are initialized with their specific constructors.
// Custom implementations of cubic.EOSType are handled by the default case, which populates
// the configuration with the substance's properties.
func (s *Substance) CubicConfig(Type cubic.EOSType, T, P, R float64) *cubic.EOSCfg {
	tc := s.Critical.Tc
	pc := s.Critical.Pc
	switch Type.(type) {
	case *cubic.VdW:
		return cubic.NewvdWCfg(T, P, tc, pc, R)
	case *cubic.RK:
		return cubic.NewRKCfg(T, P, tc, pc, R)
	case *cubic.SRK:
		return cubic.NewSRKCfg(T, P, tc, pc, s.Acentric, R)
	case *cubic.PR:
		return cubic.NewPRCfg(T, P, tc, pc, s.Acentric, R)
	default:
		return &cubic.EOSCfg{
			Type:     Type,
			T:        T,
			P:        P,
			Tc:       tc,
			Pc:       pc,
			Acentric: s.Acentric,
			R:        R,
		}
	}
}

// Vsat calculates the saturated liquid molar volume at the given temperature using the Rackett equation.
// Temperature must be in Kelvin.
func (s *Substance) Vsat(T float64) (float64, error) {
	if T <= 0 {
		return 0, zfactor.ErrTemp
	}

	tr := T / s.Critical.Tc

	return liquids.Vsat(s.Critical.Vc, s.Critical.Zc, tr)
}

// ReducedDensity calculates the reduced density (rho_r) of the substance at the given
// temperature (K) and pressure (bar) using the Lydersen chart correlation.
//
// It returns an error if the temperature is non-positive, pressure is negative,
// or if the state point is outside the range of the Lydersen chart.
func (s *Substance) ReducedDensity(T, P float64) (float64, error) {
	if T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if P < 0 {
		return 0, zfactor.ErrPressure
	}

	tr := T / s.Critical.Tc
	pr := P / s.Critical.Pc

	return liquids.ReducedDensity(tr, pr)
}

// ResidualEnthalpy calculates the dimensionless residual enthalpy H^R / (R * Tc)
// at the given temperature (K) and pressure (bar) using the Abbott (Virial) correlations.
//
// It returns an error if the temperature is non-positive or pressure is non-positive.
func (s *Substance) ResidualEnthalpy(T, P float64) (float64, error) {
	if T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if P <= 0 {
		return 0, zfactor.ErrPressure
	}
	Tr := T / s.Critical.Tc
	Pr := P / s.Critical.Pc

	return abbott.ResidualEnthalpy(Tr, Pr, s.Acentric)
}

// ResidualEntropy calculates the dimensionless residual entropy S^R / R
// at the given temperature (K) and pressure (bar) using the Abbott (Virial) correlations.
//
// It returns an error if the temperature is non-positive or pressure is non-positive.
func (s *Substance) ResidualEntropy(T, P float64) (float64, error) {
	if T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if P <= 0 {
		return 0, zfactor.ErrPressure
	}
	Tr := T / s.Critical.Tc
	Pr := P / s.Critical.Pc

	return abbott.ResidualEntropy(Tr, Pr, s.Acentric)
}
