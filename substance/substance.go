// package substance contains the characteristic properties of pure
// species
package substance

import (
	"github.com/rickykimani/zfactor"
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
// p is the pressure in bar.
// t is the temperature in Kelvin.
func (s *Substance) LeeKesler(p, t float64, property leekesler.Property) (float64, error) {
	pr := p / s.Critical.Pc
	tr := t / s.Critical.Tc

	c := leekesler.Correlation(property)

	v0, v1, err := c.At(pr, tr)
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
func (s *Substance) CubicConfig(Type cubic.EOSType, t, p, r float64) *cubic.EOSCfg {
	tc := s.Critical.Tc
	pc := s.Critical.Pc
	switch Type.(type) {
	case *cubic.VdW:
		return cubic.NewvdWCfg(t, p, tc, pc, r)
	case *cubic.RK:
		return cubic.NewRKCfg(t, p, tc, pc, r)
	case *cubic.SRK:
		return cubic.NewSRKCfg(t, p, tc, pc, s.Acentric, r)
	case *cubic.PR:
		return cubic.NewPRCfg(t, p, tc, pc, s.Acentric, r)
	default:
		return &cubic.EOSCfg{
			Type:     Type,
			T:        t,
			P:        p,
			Tc:       tc,
			Pc:       pc,
			Acentric: s.Acentric,
			R:        r,
		}
	}
}

// Vsat calculates the saturated liquid molar volume at the given temperature using the Rackett equation.
// Temperature must be in Kelvin.
func (s *Substance) Vsat(Temperature float64) (float64, error) {
	if Temperature <= 0 {
		return 0, zfactor.ErrTemp
	}

	tr := Temperature / s.Critical.Tc

	return liquids.Vsat(s.Critical.Vc, s.Critical.Zc, tr)
}

// ReducedDensity calculates the reduced density (rho_r) of the substance at the given
// temperature (K) and pressure (bar) using the Lydersen chart correlation.
//
// It returns an error if the temperature is non-positive, pressure is negative,
// or if the state point is outside the range of the Lydersen chart.
func (s *Substance) ReducedDensity(Temperature, Pressure float64) (float64, error) {
	if Temperature <= 0 {
		return 0, zfactor.ErrTemp
	}
	if Pressure < 0 {
		return 0, zfactor.ErrPressure
	}

	tr := Temperature / s.Critical.Tc
	pr := Pressure / s.Critical.Pc

	return liquids.ReducedDensity(tr, pr)
}
