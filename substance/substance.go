// package substance contains the characteristic properties of pure
// species
package substance

import (
	"github.com/rickykimani/zfactor/cubic"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
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
func (s Substance) LeeKesler(p, t float64, property leekesler.Property) (float64, error) {
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

func (s Substance) VdWCfg(t, p, r float64) *cubic.EOSCfg {
	return cubic.NewvdWCfg(t, p, s.Critical.Tc, s.Critical.Pc, r)
}
func (s Substance) RKCfg(t, p, r float64) *cubic.EOSCfg {
	return cubic.NewRKCfg(t, p, s.Critical.Tc, s.Critical.Pc, r)
}
func (s Substance) SRKCfg(t, p, r float64) *cubic.EOSCfg {
	return cubic.NewSRKCfg(t, p, s.Critical.Tc, s.Critical.Pc, s.Acentric, r)
}
func (s Substance) PRCfg(t, p, r float64) *cubic.EOSCfg {
	return cubic.NewPRCfg(t, p, s.Critical.Tc, s.Critical.Pc, s.Acentric, r)
}
