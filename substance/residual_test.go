package substance_test

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/substance"
)

// R in bar*cm^3/(mol*K), the convention the equations of state are used with.
const residualR = 10 * zfactor.RSI

// Table 13.7 of Smith, Van Ness and Abbott, reached through the substance
// methods rather than by assembling the dimensionless state by hand.
//
// The cubic package already checks these numbers against its own functions.
// Running them again through here is what confirms the wrappers pass the right
// state along: a substance whose critical properties or acentric factor went
// astray on the way in would still satisfy every test in that package.
func TestCubicResidualPropertiesAgainstTable13_7(t *testing.T) {
	const (
		bookR = 8.314 // J/(mol*K), as the book quotes it
		T     = 500.0 // K
		P     = 50.0  // bar
	)

	butane := substance.NButane
	Tc := butane.Critical.Tc

	cases := []struct {
		name  string
		eos   cubic.EOSType
		wantH float64 // J/mol
		wantS float64 // J/(mol*K)
	}{
		{"vdW", &cubic.VdW{}, -3937, -5.424},
		{"RK", &cubic.RK{}, -4505, -6.546},
		{"SRK", &cubic.SRK{}, -4824, -7.413},
		{"PR", &cubic.PR{}, -4988, -7.426},
	}

	args := zfactor.Args{T: T, P: P, R: residualR}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			hr, err := butane.CubicResidualEnthalpy(c.eos, cubic.VaporPhase, args)
			if err != nil {
				t.Fatalf("CubicResidualEnthalpy: %v", err)
			}

			sr, err := butane.CubicResidualEntropy(c.eos, cubic.VaporPhase, args)
			if err != nil {
				t.Fatalf("CubicResidualEntropy: %v", err)
			}

			gotH := hr * bookR * Tc
			gotS := sr * bookR

			// A tenth of a percent, which the rounding in the published
			// values accounts for several times over.
			const tolerance = 1e-3

			if relative := math.Abs(gotH/c.wantH - 1); relative > tolerance {
				t.Errorf("H^R: got %.6g J/mol, book gives %.6g, off by %.3g%%",
					gotH, c.wantH, 100*relative)
			}

			if relative := math.Abs(gotS/c.wantS - 1); relative > tolerance {
				t.Errorf("S^R: got %.6g J/(mol K), book gives %.6g, off by %.3g%%",
					gotS, c.wantS, 100*relative)
			}
		})
	}
}

// The residual Gibbs energy ties the other two together,
//
//	ln phi = H^R/(R T) - S^R/R
//
// so the three methods have to agree with each other when called separately.
// This checks that they are all describing the same state, which is the thing
// that could go wrong in a wrapper even when the underlying functions are
// correct.
func TestCubicDepartureMethodsDescribeTheSameState(t *testing.T) {
	butane := substance.NButane
	Tc := butane.Critical.Tc

	states := []struct {
		name string
		T, P float64
	}{
		{"gas", 500, 50},
		{"low pressure", 500, 1},
		{"near critical", 430, 40},
	}

	equations := map[string]cubic.EOSType{
		"vdW": &cubic.VdW{},
		"RK":  &cubic.RK{},
		"SRK": &cubic.SRK{},
		"PR":  &cubic.PR{},
	}

	for name, eos := range equations {
		for _, s := range states {
			t.Run(name+" "+s.name, func(t *testing.T) {
				args := zfactor.Args{T: s.T, P: s.P, R: residualR}

				hr, err := butane.CubicResidualEnthalpy(eos, cubic.VaporPhase, args)
				if err != nil {
					t.Fatalf("CubicResidualEnthalpy: %v", err)
				}

				sr, err := butane.CubicResidualEntropy(eos, cubic.VaporPhase, args)
				if err != nil {
					t.Fatalf("CubicResidualEntropy: %v", err)
				}

				lnPhi, err := butane.CubicLogFugacity(eos, cubic.VaporPhase, args)
				if err != nil {
					t.Fatalf("CubicLogFugacity: %v", err)
				}

				// The enthalpy is normalised by R Tc and the identity is in
				// terms of H^R/(R T).
				overRT := hr / (s.T / Tc)

				const tolerance = 1e-12

				if diff := math.Abs((overRT - sr) - lnPhi); diff > tolerance {
					t.Errorf("H^R/RT - S^R/R = %.15g, ln phi = %.15g, differ by %.3g",
						overRT-sr, lnPhi, diff)
				}
			})
		}
	}
}

// Inside the two-phase region the equation has three roots, and the liquid and
// the vapour have genuinely different residual properties. Naming the phase has
// to select between them, or the argument is doing nothing.
func TestCubicResidualPropertiesDifferByPhase(t *testing.T) {
	butane := substance.NButane

	// Below the critical temperature and near the saturation pressure, so all
	// three roots are real.
	T := 0.85 * butane.Critical.Tc

	eos := &cubic.PR{}

	cfg := butane.CubicConfig(eos, zfactor.Args{T: T, R: residualR})

	pSat, err := cubic.SaturationPressure(cfg, T)
	if err != nil {
		t.Fatalf("SaturationPressure at T=%g: %v", T, err)
	}

	args := zfactor.Args{T: T, P: pSat, R: residualR}

	// The configuration above carries no pressure, since SaturationPressure
	// is given the temperature and finds the pressure itself. Solving for a
	// phase needs both.
	satCfg := butane.CubicConfig(eos, args)

	// Confirm the premise: the state really does have two phases to choose
	// between, or this would be comparing a root against itself.
	vaporState, err := cubic.SolvePhase(satCfg, cubic.VaporPhase)
	if err != nil {
		t.Fatalf("SolvePhase vapor: %v", err)
	}

	liquidState, err := cubic.SolvePhase(satCfg, cubic.LiquidPhase)
	if err != nil {
		t.Fatalf("SolvePhase liquid: %v", err)
	}

	if vaporState.Z == liquidState.Z {
		t.Fatalf("both phases gave Z=%g, so the state has a single root and this test is vacuous",
			vaporState.Z)
	}

	vaporH, err := butane.CubicResidualEnthalpy(eos, cubic.VaporPhase, args)
	if err != nil {
		t.Fatalf("vapor enthalpy: %v", err)
	}

	liquidH, err := butane.CubicResidualEnthalpy(eos, cubic.LiquidPhase, args)
	if err != nil {
		t.Fatalf("liquid enthalpy: %v", err)
	}

	// Condensing releases energy, so the liquid is the more negative of the
	// two by the enthalpy of vaporisation.
	if liquidH >= vaporH {
		t.Errorf("liquid H^R/(R Tc) = %.6g is not below the vapour %.6g", liquidH, vaporH)
	}

	vaporS, err := butane.CubicResidualEntropy(eos, cubic.VaporPhase, args)
	if err != nil {
		t.Fatalf("vapor entropy: %v", err)
	}

	liquidS, err := butane.CubicResidualEntropy(eos, cubic.LiquidPhase, args)
	if err != nil {
		t.Fatalf("liquid entropy: %v", err)
	}

	// The liquid is the more ordered phase, so its residual entropy is the
	// more negative.
	if liquidS >= vaporS {
		t.Errorf("liquid S^R/R = %.6g is not below the vapour %.6g", liquidS, vaporS)
	}

	// At the saturation pressure the two phases are in equilibrium, which is
	// the condition SaturationPressure solves: equal fugacity.
	vaporPhi, err := butane.CubicLogFugacity(eos, cubic.VaporPhase, args)
	if err != nil {
		t.Fatalf("vapor fugacity: %v", err)
	}

	liquidPhi, err := butane.CubicLogFugacity(eos, cubic.LiquidPhase, args)
	if err != nil {
		t.Fatalf("liquid fugacity: %v", err)
	}

	if diff := math.Abs(vaporPhi - liquidPhi); diff > 1e-6 {
		t.Errorf("fugacities at the saturation pressure differ by %.3g: vapor %.8g, liquid %.8g",
			diff, vaporPhi, liquidPhi)
	}
}

func TestCubicResidualPropertiesRejectBadArgs(t *testing.T) {
	butane := substance.NButane
	eos := &cubic.PR{}

	cases := []struct {
		name string
		args zfactor.Args
		want error
	}{
		{"no temperature", zfactor.Args{T: 0, P: 50, R: residualR}, zfactor.ErrTemp},
		{"negative temperature", zfactor.Args{T: -1, P: 50, R: residualR}, zfactor.ErrTemp},
		{"no pressure", zfactor.Args{T: 500, P: 0, R: residualR}, zfactor.ErrPressure},
		{"no gas constant", zfactor.Args{T: 500, P: 50, R: 0}, zfactor.ErrUniversalConst},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := butane.CubicResidualEnthalpy(eos, cubic.VaporPhase, c.args); !errors.Is(err, c.want) {
				t.Errorf("CubicResidualEnthalpy: got %v, want %v", err, c.want)
			}
			if _, err := butane.CubicResidualEntropy(eos, cubic.VaporPhase, c.args); !errors.Is(err, c.want) {
				t.Errorf("CubicResidualEntropy: got %v, want %v", err, c.want)
			}
			if _, err := butane.CubicLogFugacity(eos, cubic.VaporPhase, c.args); !errors.Is(err, c.want) {
				t.Errorf("CubicLogFugacity: got %v, want %v", err, c.want)
			}
		})
	}
}
