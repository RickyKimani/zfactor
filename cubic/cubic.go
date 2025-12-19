package cubic

import (
	"fmt"
	"math"
	"slices"

	"github.com/rickykimani/zfactor"
)


// Params represents the substance agnostic variables in any
// cubic equation of state
type Params struct {
	Sigma   float64 //σ
	Epsilon float64 //ε
	Omega   float64 //Ω
	Psi     float64 //Ψ
}

// EOSType defines what makes up an equation of state
type EOSType interface {
	Alpha(tr, w float64) float64 //α(Tr, ω)
	Params() *Params
}

// VolumeResult contains the results of solving the cubic equation of state for volume.
type VolumeResult struct {
	A       float64       // The a(T) parameter value
	B       float64       // The b parameter value
	Volumes [3]complex128 // The roots of the cubic equation (molar volumes)
}

// Clean returns the real roots of the volume equation, sorted in ascending order.
// The smallest root corresponds to the liquid phase, and the largest to the vapor phase.
func (vr *VolumeResult) Clean() []float64 {
	res := make([]float64, 0, 3)
	for _, value := range vr.Volumes {
		if math.Abs(imag(value)) < 1e-9 {
			res = append(res, real(value))
		}
	}
	slices.Sort(res)
	return res
}

// String implements fmt.Stringer for VolumeResult.
func (vr *VolumeResult) String() string {
	return fmt.Sprintf("VolumeResult{A: %g, B: %g, Volumes: %v}", vr.A, vr.B, vr.Volumes)
}

// PressureResult contains the calculated pressure and intermediate parameters.
type PressureResult struct {
	A float64 // The a(T) parameter value
	B float64 // The b parameter value
	P float64 // The calculated pressure
}

// String implements fmt.Stringer for PressureResult.
func (pr *PressureResult) String() string {
	return fmt.Sprintf("PressureResult{A: %g, B: %g, P: %g}", pr.A, pr.B, pr.P)
}

// EOSCfg holds the configuration and state variables for an Equation of State calculation.
type EOSCfg struct {
	Type     EOSType // The type of cubic equation of state (e.g., VdW, RK, SRK, PR)
	T        float64 // Absolute temperature
	P        float64 // Pressure
	Tc       float64 // Critical temperature
	Pc       float64 // Critical pressure
	Acentric float64 // Acentric factor (ω) - dimensionless
	R        float64 // Universal gas constant in consistent units
}

// calculateb calculates the b parameter
func calculateb(omega, r, tc, pc float64) float64 {
	return omega * r * tc / pc
}

// calculatea calculates the a(T) parameter
func calculatea(psi, alpha, r, tc, pc float64) float64 {
	return psi * alpha * r * r * tc * tc / pc
}

// SolveForVolume solves the cubic equation of state for molar volume given the configuration.
// It returns the calculated parameters a and b, and the three roots of the cubic equation.
// Returns an error if input parameters are invalid (e.g. non-positive temperature).
func SolveForVolume(cfg *EOSCfg) (*VolumeResult, error) {
	if cfg.T <= 0 {
		return nil, zfactor.TempErr
	}
	if cfg.P <= 0 {
		return nil, zfactor.PressErr
	}

	if cfg.Pc <= 0 || cfg.Tc <= 0 {
		return nil, zfactor.CriticalPropErr
	}

	if cfg.R <= 0 {
		return nil, zfactor.UniversalConstErr
	}

	tr := cfg.T / cfg.Tc

	alpha := cfg.Type.Alpha(tr, cfg.Acentric)

	sigma := cfg.Type.Params().Sigma
	epsilon := cfg.Type.Params().Epsilon
	omega := cfg.Type.Params().Omega
	psi := cfg.Type.Params().Psi

	a := calculatea(psi, alpha, cfg.R, cfg.Tc, cfg.Pc)
	b := calculateb(omega, cfg.R, cfg.Tc, cfg.Pc)

	//eV^3 + fV^2 + gV + h = 0
	x := epsilon + sigma
	y := epsilon * sigma
	v_ig := cfg.R * cfg.Tc / cfg.Pc

	e := 1.0
	f := b*(x-1) - v_ig
	g := b*((y-x)*b-(x*v_ig)) + a/cfg.P
	h := -y*b*b*(b+v_ig) - a*b/cfg.P

	solution, err := zfactor.SolveCubic(e, f, g, h)
	if err != nil {
		return nil, fmt.Errorf("failed to solve cubic: %w", err)
	}

	return &VolumeResult{
		A:       a,
		B:       b,
		Volumes: solution,
	}, nil

}

// Pressure calculates the pressure for a given molar volume and configuration.
// It returns the calculated pressure and parameters a and b.
// Returns an error if input parameters are invalid.
func Pressure(cfg *EOSCfg, volume float64) (*PressureResult, error) {
	if cfg.T <= 0 {
		return nil, zfactor.TempErr
	}
	if cfg.Pc <= 0 || cfg.Tc <= 0 {
		return nil, zfactor.CriticalPropErr
	}

	if cfg.R <= 0 {
		return nil, zfactor.UniversalConstErr
	}
	tr := cfg.T / cfg.Tc

	alpha := cfg.Type.Alpha(tr, cfg.Acentric)

	sigma := cfg.Type.Params().Sigma
	epsilon := cfg.Type.Params().Epsilon
	omega := cfg.Type.Params().Omega
	psi := cfg.Type.Params().Psi

	a := calculatea(psi, alpha, cfg.R, cfg.Tc, cfg.Pc)
	b := calculateb(omega, cfg.R, cfg.Tc, cfg.Pc)
	v := volume

	first := cfg.R * cfg.T / (v - b)
	second := a / ((v + epsilon*b) * (v + sigma*b))

	p := first - second

	return &PressureResult{
		A: a,
		B: b,
		P: p,
	}, nil
}

type vdW struct{}

func (*vdW) Alpha(tr, w float64) float64 {
	return 1.0
}

func (*vdW) Params() *Params {
	return &Params{
		Sigma:   0,
		Epsilon: 0,
		Omega:   1.0 / 8.0,
		Psi:     27.0 / 64.0,
	}
}

// NewvdWCfg creates a configuration for the van der Waals cubic equation of state
func NewvdWCfg(T, P, Tc, Pc, R float64) *EOSCfg {
	return &EOSCfg{
		Type:     &vdW{},
		T:        T,
		P:        P,
		Tc:       Tc,
		Pc:       Pc,
		Acentric: 0,
		R:        R,
	}
}

type rk struct{}

func (*rk) Alpha(tr, w float64) float64 {
	return 1 / math.Sqrt(tr)
}

func (*rk) Params() *Params {
	return &Params{
		Sigma:   1,
		Epsilon: 0,
		Omega:   0.08664,
		Psi:     0.42728,
	}
}

// NewRKCfg creates a configuration for the Redlich-Kwong cubic equation of state
func NewRKCfg(T, P, Tc, Pc, R float64) *EOSCfg {
	return &EOSCfg{
		Type:     &rk{},
		T:        T,
		P:        P,
		Tc:       Tc,
		Pc:       Pc,
		Acentric: 0,
		R:        R,
	}
}

type srk struct{}

func (*srk) Alpha(tr, w float64) float64 {
	a := 0.480 + 1.574*w - 0.716*w*w
	b := 1 - math.Sqrt(tr)
	c := 1 + a*b
	return c * c
}

func (*srk) Params() *Params {
	return &Params{
		Sigma:   1,
		Epsilon: 0,
		Omega:   0.08664,
		Psi:     0.42728,
	}
}

// NewSRKCfg creates a configuration for the Soave-Redlich-Kwong cubic equation of state
func NewSRKCfg(T, P, Tc, Pc, W, R float64) *EOSCfg {
	return &EOSCfg{
		Type:     &srk{},
		T:        T,
		P:        P,
		Tc:       Tc,
		Pc:       Pc,
		Acentric: W,
		R:        R,
	}
}

type pr struct{}

func (*pr) Alpha(tr, w float64) float64 {
	a := 0.37464 + 1.54226*w - 0.26992*w*w
	b := 1 - math.Sqrt(tr)
	c := 1 + a*b
	return c * c
}

func (*pr) Params() *Params {
	return &Params{
		Sigma:   1 + math.Sqrt2,
		Epsilon: 1 - math.Sqrt2,
		Omega:   0.07780,
		Psi:     0.45724,
	}
}

// NewPRCfg creates a configuration for the Peng-Robinson cubic equation of state
func NewPRCfg(T, P, Tc, Pc, W, R float64) *EOSCfg {
	return &EOSCfg{
		Type:     &pr{},
		T:        T,
		P:        P,
		Tc:       Tc,
		Pc:       Pc,
		Acentric: W,
		R:        R,
	}
}
