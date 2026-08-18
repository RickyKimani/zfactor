package cubic_test

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/substance"
)

// A state below the critical temperature at its saturation pressure has three
// real roots, and the two phases must select the outer ones. This is the case
// the phase argument exists for.
func TestSolvePhaseSelectsTheOuterRoots(t *testing.T) {
	sub := substance.NButane
	T := 0.85 * sub.Critical.Tc

	base := cubic.NewPRCfg(T, 0, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	pSat, err := cubic.SaturationPressure(base, T)
	if err != nil {
		t.Fatalf("SaturationPressure: %v", err)
	}

	cfg := *base
	cfg.P = pSat

	// The roots this is selecting among, so the comparison below is against
	// the volumes themselves rather than against another belief about them.
	result, err := cubic.SolveForVolume(&cfg)
	if err != nil {
		t.Fatalf("SolveForVolume: %v", err)
	}

	roots := result.Clean()
	if len(roots) != 3 {
		t.Fatalf("state has %d real roots, want 3; this test needs a two-phase state", len(roots))
	}

	RT := cfg.R * T

	vapor, err := cubic.SolvePhase(&cfg, cubic.VaporPhase)
	if err != nil {
		t.Fatalf("SolvePhase vapor: %v", err)
	}

	liquid, err := cubic.SolvePhase(&cfg, cubic.LiquidPhase)
	if err != nil {
		t.Fatalf("SolvePhase liquid: %v", err)
	}

	wantVapor := pSat * roots[2] / RT
	wantLiquid := pSat * roots[0] / RT

	if math.Abs(vapor.Z-wantVapor) > 1e-12 {
		t.Errorf("vapor Z = %.12g, want %.12g from the largest root", vapor.Z, wantVapor)
	}

	if math.Abs(liquid.Z-wantLiquid) > 1e-12 {
		t.Errorf("liquid Z = %.12g, want %.12g from the smallest root", liquid.Z, wantLiquid)
	}

	// The middle root is a state of negative compressibility with respect to
	// volume and is not physical, so neither phase should land on it.
	middle := pSat * roots[1] / RT

	if math.Abs(vapor.Z-middle) < 1e-12 || math.Abs(liquid.Z-middle) < 1e-12 {
		t.Error("a phase selected the middle root")
	}

	// A and B describe the equation at the state, not the phase, so they are
	// the same for both.
	if vapor.A != liquid.A || vapor.B != liquid.B {
		t.Errorf("A and B differ between phases: vapor (%g, %g), liquid (%g, %g)",
			vapor.A, vapor.B, liquid.A, liquid.B)
	}
}

// Where the equation has a single real root there is only one way for the
// fluid to be, and both phases name it.
func TestSolvePhaseWithASingleRootGivesTheSameForBoth(t *testing.T) {
	sub := substance.NButane

	// Well above the critical temperature at a low pressure, so there is one
	// real root.
	cfg := cubic.NewPRCfg(700, 1, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	result, err := cubic.SolveForVolume(cfg)
	if err != nil {
		t.Fatalf("SolveForVolume: %v", err)
	}

	if roots := result.Clean(); len(roots) != 1 {
		t.Fatalf("state has %d real roots, want 1", len(roots))
	}

	vapor, err := cubic.SolvePhase(cfg, cubic.VaporPhase)
	if err != nil {
		t.Fatalf("SolvePhase vapor: %v", err)
	}

	liquid, err := cubic.SolvePhase(cfg, cubic.LiquidPhase)
	if err != nil {
		t.Fatalf("SolvePhase liquid: %v", err)
	}

	if vapor.Z != liquid.Z {
		t.Errorf("single root gave Z=%g for vapor and Z=%g for liquid", vapor.Z, liquid.Z)
	}
}

// The dimensionless state is the triple the property functions take, so it
// should match what forming it by hand gives.
func TestSolvePhaseReducesToTheDimensionlessGroups(t *testing.T) {
	sub := substance.NButane

	const T, P = 500.0, 50.0

	cfg := cubic.NewRKCfg(T, P, sub.Critical.Tc, sub.Critical.Pc, residualR)

	result, err := cubic.SolveForVolume(cfg)
	if err != nil {
		t.Fatalf("SolveForVolume: %v", err)
	}

	state, err := cubic.SolvePhase(cfg, cubic.VaporPhase)
	if err != nil {
		t.Fatalf("SolvePhase: %v", err)
	}

	RT := residualR * T

	wantA := result.A * P / (RT * RT)
	wantB := result.B * P / RT

	if math.Abs(state.A-wantA) > 1e-15 {
		t.Errorf("A = %.15g, want %.15g", state.A, wantA)
	}

	if math.Abs(state.B-wantB) > 1e-15 {
		t.Errorf("B = %.15g, want %.15g", state.B, wantB)
	}
}

func TestSolvePhaseRejectsBadConfig(t *testing.T) {
	sub := substance.NButane
	valid := cubic.NewPRCfg(500, 50, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	noTemperature := *valid
	noTemperature.T = 0

	noPressure := *valid
	noPressure.P = 0

	noGasConstant := *valid
	noGasConstant.R = 0

	cases := []struct {
		name string
		cfg  *cubic.EOSCfg
	}{
		{"nil config", nil},
		{"no temperature", &noTemperature},
		{"no pressure", &noPressure},
		{"no gas constant", &noGasConstant},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := cubic.SolvePhase(c.cfg, cubic.VaporPhase); err == nil {
				t.Error("SolvePhase returned no error")
			}
		})
	}
}

func TestSolvePhaseRejectsUnknownPhase(t *testing.T) {
	sub := substance.NButane
	cfg := cubic.NewPRCfg(500, 50, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	if _, err := cubic.SolvePhase(cfg, cubic.Phase(99)); err == nil {
		t.Error("SolvePhase accepted an unknown phase")
	}
}

// A phase is named in error messages, so the names have to be right. The
// unknown case is included because a Phase is an integer and a caller can
// produce one that is not a phase at all.
func TestPhaseString(t *testing.T) {
	cases := []struct {
		phase cubic.Phase
		want  string
	}{
		{cubic.StablePhase, "stable"},
		{cubic.VaporPhase, "vapor"},
		{cubic.LiquidPhase, "liquid"},
		{cubic.Phase(7), "Phase(7)"},
	}

	for _, c := range cases {
		if got := c.phase.String(); got != c.want {
			t.Errorf("Phase(%d).String() = %q, want %q", int(c.phase), got, c.want)
		}
	}
}

// ErrNoRealRoot is documented as the sentinel for a state the equation cannot
// place, so it has to be reachable through errors.Is for a caller to act on.
func TestErrNoRealRootIsASentinel(t *testing.T) {
	if !errors.Is(cubic.ErrNoRealRoot, cubic.ErrNoRealRoot) {
		t.Error("ErrNoRealRoot does not match itself")
	}

	if cubic.ErrNoRealRoot.Error() == "" {
		t.Error("ErrNoRealRoot has no message")
	}
}

// Stability can be characterised two ways, and they must agree: the stable
// phase is the one of lower fugacity, and it is also the liquid above the
// saturation pressure and the vapour below it.
//
// SolvePhase uses the first because it needs no iteration. This checks it
// against the second, which is the rule the diagram code used before it was
// replaced, so the two independent characterisations are held against each
// other over states that genuinely have a phase to choose.
func TestStablePhaseAgreesWithTheSaturationPressureRule(t *testing.T) {
	sub := substance.NButane
	Tc := sub.Critical.Tc

	// Subcritical, where three roots are possible.
	for _, reduced := range []float64{0.7, 0.8, 0.9, 0.95} {
		T := reduced * Tc

		base := cubic.NewPRCfg(T, 0, Tc, sub.Critical.Pc, sub.Acentric, residualR)

		pSat, err := cubic.SaturationPressure(base, T)
		if err != nil {
			t.Fatalf("SaturationPressure at Tr=%g: %v", reduced, err)
		}

		// Straddle the saturation pressure. Far from it the equation has a
		// single root and the choice is moot, which is itself worth covering.
		for _, fraction := range []float64{0.5, 0.9, 0.99, 1.01, 1.1, 1.5} {
			cfg := *base
			cfg.P = fraction * pSat

			result, err := cubic.SolveForVolume(&cfg)
			if err != nil {
				continue
			}

			roots := result.Clean()
			if len(roots) == 0 {
				continue
			}

			state, err := cubic.SolvePhase(&cfg, cubic.StablePhase)
			if err != nil {
				t.Fatalf("SolvePhase at Tr=%g P/Psat=%g: %v", reduced, fraction, err)
			}

			// Above the saturation pressure the liquid is stable, below it
			// the vapour is.
			want := roots[len(roots)-1]
			wantName := "vapor"

			if fraction > 1 {
				want = roots[0]
				wantName = "liquid"
			}

			if state.V != want {
				t.Errorf("Tr=%g P/Psat=%g (%d roots): stable phase gave V=%.6g, the saturation rule says %s at V=%.6g",
					reduced, fraction, len(roots), state.V, wantName, want)
			}
		}
	}
}

// At the saturation pressure the two phases are in equilibrium and their
// fugacities are equal, so the comparison is a genuine tie. The documented
// convention is the vapour, rather than whichever side the floating-point
// comparison falls on.
func TestStablePhaseBreaksTheSaturationTieTowardVapor(t *testing.T) {
	sub := substance.NButane
	T := 0.85 * sub.Critical.Tc

	base := cubic.NewPRCfg(T, 0, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	pSat, err := cubic.SaturationPressure(base, T)
	if err != nil {
		t.Fatalf("SaturationPressure: %v", err)
	}

	cfg := *base
	cfg.P = pSat

	result, err := cubic.SolveForVolume(&cfg)
	if err != nil {
		t.Fatalf("SolveForVolume: %v", err)
	}

	roots := result.Clean()
	if len(roots) != 3 {
		t.Fatalf("state at the saturation pressure has %d real roots, want 3", len(roots))
	}

	state, err := cubic.SolvePhase(&cfg, cubic.StablePhase)
	if err != nil {
		t.Fatalf("SolvePhase: %v", err)
	}

	if vapor := roots[2]; state.V != vapor {
		t.Errorf("at the saturation pressure the tie gave V=%.6g, want the vapour at V=%.6g",
			state.V, vapor)
	}
}

// The zero value is the stable phase, so a caller who has not thought about
// phase gets the physical answer instead of an arbitrary root.
func TestZeroPhaseIsStable(t *testing.T) {
	var zero cubic.Phase

	if zero != cubic.StablePhase {
		t.Errorf("the zero Phase is %v, want %v", zero, cubic.StablePhase)
	}
}

// The volume carried alongside the dimensionless groups is the root itself, so
// reducing it should give back the compressibility factor.
func TestPhaseStateVolumeAndCompressibilityAgree(t *testing.T) {
	sub := substance.NButane

	const T, P = 500.0, 50.0

	cfg := cubic.NewPRCfg(T, P, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	state, err := cubic.SolvePhase(cfg, cubic.StablePhase)
	if err != nil {
		t.Fatalf("SolvePhase: %v", err)
	}

	want := P * state.V / (residualR * T)

	if math.Abs(state.Z-want) > 1e-15 {
		t.Errorf("Z = %.15g, but PV/RT from the reported volume is %.15g", state.Z, want)
	}
}
