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

	return property.Eval(tr, pr, s.Acentric)
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
	return &cubic.EOSCfg{
		Type:     Type,
		T:        args.T,
		P:        args.P,
		Tc:       s.Critical.Tc,
		Pc:       s.Critical.Pc,
		Acentric: s.Acentric,
		R:        args.R,
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

// CubicResidualEnthalpy calculates the dimensionless residual enthalpy
// H^R / (R * Tc) at the given temperature (K) and pressure (bar) from a cubic
// equation of state.
//
// This is the same quantity AbbottResidualEnthalpy and the Lee-Kesler tables
// return, by a different route: the equation of state is solved at the state
// and the departure functions evaluated from its own parameters, rather than
// read from a generalized correlation.
//
// A state may admit a liquid root and a vapour root, and it has a residual
// enthalpy for each, so the phase is named rather than assumed. Where only one
// real root exists both phases return it.
//
// Required Args:
//   - T: Temperature in Kelvin
//   - P: Pressure in bar
//   - R: Gas constant, in units consistent with the critical pressure
//
// It returns an error if the temperature or pressure is non-positive, or if the
// equation has no real root at the state.
func (s *Substance) CubicResidualEnthalpy(Type cubic.EOSType, phase cubic.Phase, args zfactor.Args) (float64, error) {
	cfg, state, err := s.cubicPhase(Type, phase, args)
	if err != nil {
		return 0, err
	}

	return cubic.ResidualEnthalpy(cfg, state.Z, state.A, state.B)
}

// CubicResidualEntropy calculates the dimensionless residual entropy S^R / R at
// the given temperature (K) and pressure (bar) from a cubic equation of state.
//
// The arguments carry the same meaning as in CubicResidualEnthalpy.
func (s *Substance) CubicResidualEntropy(Type cubic.EOSType, phase cubic.Phase, args zfactor.Args) (float64, error) {
	cfg, state, err := s.cubicPhase(Type, phase, args)
	if err != nil {
		return 0, err
	}

	return cubic.ResidualEntropy(cfg, state.Z, state.A, state.B)
}

// CubicLogFugacity calculates the natural logarithm of the fugacity
// coefficient at the given temperature (K) and pressure (bar) from a cubic
// equation of state.
//
// It completes the set: the three departure properties a cubic equation
// supplies are the residual enthalpy, the residual entropy and this, and
// ln phi is the residual Gibbs energy G^R/(R T) that ties the other two
// together.
//
// The arguments carry the same meaning as in CubicResidualEnthalpy.
func (s *Substance) CubicLogFugacity(Type cubic.EOSType, phase cubic.Phase, args zfactor.Args) (float64, error) {
	cfg, state, err := s.cubicPhase(Type, phase, args)
	if err != nil {
		return 0, err
	}

	return cubic.LogFugacity(cfg, state.Z, state.A, state.B), nil
}

// cubicPhase builds the configuration for a cubic equation of state at a state
// and solves it for one phase, which is the work every method above shares.
func (s *Substance) cubicPhase(
	Type cubic.EOSType,
	phase cubic.Phase,
	args zfactor.Args,
) (*cubic.EOSCfg, *cubic.PhaseState, error) {
	if args.T <= 0 {
		return nil, nil, zfactor.ErrTemp
	}
	if args.P <= 0 {
		return nil, nil, zfactor.ErrPressure
	}
	if args.R <= 0 {
		return nil, nil, zfactor.ErrUniversalConst
	}

	cfg := s.CubicConfig(Type, args)

	state, err := cubic.SolvePhase(cfg, phase)
	if err != nil {
		return nil, nil, err
	}

	return cfg, state, nil
}
