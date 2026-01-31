// package substance contains the characteristic properties of pure
// species
package substance

import (
	"fmt"

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
//
// Required Args:
//   - T: Temperature in Kelvin
//   - P: Pressure in bar
func (s *Substance) LeeKesler(args zfactor.Args, property leekesler.Property) (float64, error) {
	pr := args.P / s.Critical.Pc
	tr := args.T / s.Critical.Tc

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
//
// Required Args:
//   - T: Temperature
//   - P: Pressure
//   - R: Gas Constant
func (s *Substance) CubicConfig(Type cubic.EOSType, args zfactor.Args) *cubic.EOSCfg {
	tc := s.Critical.Tc
	pc := s.Critical.Pc
	switch Type.(type) {
	case *cubic.VdW:
		return cubic.NewvdWCfg(args.T, args.P, tc, pc, args.R)
	case *cubic.RK:
		return cubic.NewRKCfg(args.T, args.P, tc, pc, args.R)
	case *cubic.SRK:
		return cubic.NewSRKCfg(args.T, args.P, tc, pc, s.Acentric, args.R)
	case *cubic.PR:
		return cubic.NewPRCfg(args.T, args.P, tc, pc, s.Acentric, args.R)
	default:
		return &cubic.EOSCfg{
			Type:     Type,
			T:        args.T,
			P:        args.P,
			Tc:       tc,
			Pc:       pc,
			Acentric: s.Acentric,
			R:        args.R,
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
//
// Required Args:
//   - T: Temperature in Kelvin
//   - P: Pressure in bar
func (s *Substance) ReducedDensity(args zfactor.Args) (float64, error) {
	if args.T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if args.P < 0 {
		return 0, zfactor.ErrPressure
	}

	tr := args.T / s.Critical.Tc
	pr := args.P / s.Critical.Pc

	return liquids.ReducedDensity(tr, pr)
}

// AbbottResidualEnthalpy calculates the dimensionless residual enthalpy H^R / (R * Tc)
// at the given temperature (K) and pressure (bar) using the Abbott (Virial) correlations.
//
// Required Args:
//   - T: Temperature in Kelvin
//   - P: Pressure in bar
//
// It returns an error if the temperature is non-positive or pressure is non-positive.
func (s *Substance) AbbottResidualEnthalpy(args zfactor.Args) (float64, error) {
	if args.T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if args.P <= 0 {
		return 0, zfactor.ErrPressure
	}
	Tr := args.T / s.Critical.Tc
	Pr := args.P / s.Critical.Pc

	return abbott.ResidualEnthalpy(Tr, Pr, s.Acentric)
}

// AbbottResidualEntropy calculates the dimensionless residual entropy S^R / R
// at the given temperature (K) and pressure (bar) using the Abbott (Virial) correlations.
//
// Required Args:
//   - T: Temperature in Kelvin
//   - P: Pressure in bar
//
// It returns an error if the temperature is non-positive or pressure is non-positive.
func (s *Substance) AbbottResidualEntropy(args zfactor.Args) (float64, error) {
	if args.T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if args.P <= 0 {
		return 0, zfactor.ErrPressure
	}
	Tr := args.T / s.Critical.Tc
	Pr := args.P / s.Critical.Pc

	return abbott.ResidualEntropy(Tr, Pr, s.Acentric)
}

// LeeKeslerAcentric estimates the acentric factor using the Lee-Kesler correlation.
// Use this if the substance has no defined acentric factor but has a known Normal Boiling Point (Tn).
func (s *Substance) LeeKeslerAcentric() (float64, error) {
	if s.Tn == 0 {
		return 0, fmt.Errorf("%s has no defined normal boiling point", s.Name)
	}
	return leekesler.EstimateAcentricFactor(s.Tn, s.Critical.Tc, s.Critical.Pc)
}

// LeeKeslerVaporPressure estimates the saturation vapor pressure (Psat) in bar at temperature T (K).
// It uses the Lee-Kesler correlation which internally estimates the acentric factor based on Tn.
func (s *Substance) LeeKeslerVaporPressure(T float64) (float64, error) {
	if s.Tn == 0 {
		return 0, fmt.Errorf("%s has no defined normal boiling point", s.Name)
	}
	return leekesler.VaporPressure(T, s.Tn, s.Critical.Tc, s.Critical.Pc)
}
