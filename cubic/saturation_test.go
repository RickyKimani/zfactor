package cubic_test

import (
	"fmt"
	"math"
	"strings"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/substance"
)

// fugacityAtSaturation returns the logarithms of the fugacity
// coefficients of the coexisting liquid and vapor at the given state.
//
// It reproduces the dimensionless grouping SaturationPressure works in:
// A = aP/(RT)^2 and B = bP/RT, with Z = PV/RT for each root.
func fugacityAtSaturation(
	t *testing.T,
	cfg *cubic.EOSCfg,
	temperature, pressure float64,
) (liquid, vapor float64, roots []float64) {
	t.Helper()

	state := *cfg
	state.T = temperature
	state.P = pressure

	res, err := cubic.SolveForVolume(&state)
	if err != nil {
		t.Fatalf("SolveForVolume returned an unexpected error: %v", err)
	}

	roots = res.Clean()
	if len(roots) != 3 {
		return 0, 0, roots
	}

	RT := cfg.R * temperature
	a := res.A * pressure / (RT * RT)
	b := res.B * pressure / RT

	zLiquid := pressure * roots[0] / RT
	zVapor := pressure * roots[len(roots)-1] / RT

	return cubic.LogFugacity(&state, zLiquid, a, b),
		cubic.LogFugacity(&state, zVapor, a, b),
		roots
}

// TestSaturationPressureEqualFugacity checks the condition the solver
// exists to satisfy.
//
// Two phases coexist when the fugacity of each is equal in both, so at
// the returned pressure the liquid and vapor roots of the equation of
// state must share a fugacity coefficient. This is the definition of the
// answer rather than a value read from anywhere, so it holds for every
// equation of state and every substance.
//
// A correct pressure also implies three real roots: the two-phase region
// is exactly where the isotherm is not monotonic.
func TestSaturationPressureEqualFugacity(t *testing.T) {
	const tol = 1e-7

	equations := map[string]cubic.EOSType{
		"vdW": &cubic.VdW{},
		"RK":  &cubic.RK{},
		"SRK": &cubic.SRK{},
		"PR":  &cubic.PR{},
	}

	species := map[string]*substance.Substance{
		"n-butane":       substance.NButane,
		"carbon dioxide": substance.CarbonDioxide,
		"propane":        substance.Propane,
	}

	for eosName, eos := range equations {
		for speciesName, s := range species {
			for _, tr := range []float64{0.7, 0.8, 0.9} {
				temperature := tr * s.Critical.Tc

				t.Run(fmt.Sprintf("%s/%s/Tr=%.1f", eosName, speciesName, tr), func(t *testing.T) {
					cfg := s.CubicConfig(eos, zfactor.Args{T: temperature, P: s.Critical.Pc, R: R})

					psat, err := cubic.SaturationPressure(cfg, temperature)
					if err != nil {
						t.Fatalf("SaturationPressure at Tr = %.1f returned an unexpected error: %v", tr, err)
					}

					if psat <= 0 || psat > s.Critical.Pc {
						t.Fatalf("saturation pressure %.4f bar at Tr = %.1f is outside (0, Pc]", psat, tr)
					}

					liquid, vapor, roots := fugacityAtSaturation(t, cfg, temperature, psat)

					if len(roots) != 3 {
						t.Fatalf(
							"at Tr = %.1f and the reported saturation pressure there are %d real roots; want 3",
							tr, len(roots),
						)
					}

					if math.Abs(liquid-vapor) > tol {
						t.Errorf(
							"at Tr = %.1f: ln(phi) is %.10f for the liquid and %.10f for the vapor",
							tr, liquid, vapor,
						)
					}
				})
			}
		}
	}
}

// TestSaturationPressureAgainstExperiment compares the predicted vapor
// pressure of n-butane at 350 K with the measured value quoted in
// Example 3.9 of Smith, Van Ness & Abbott.
//
// The equations are ordered by how well they reproduce it. That ordering
// is the point of the comparison: van der Waals is a qualitative model,
// Redlich/Kwong improves on it, and the two temperature-dependent forms
// were fitted with vapor pressure in mind, so Peng/Robinson lands within
// a fraction of a percent. A change that disturbed the ordering would
// indicate a defect in one of the alpha functions.
func TestSaturationPressureAgainstExperiment(t *testing.T) {
	const (
		T          = 350.0
		measured   = 9.4573 // bar
		acceptable = 0.02   // for the temperature-dependent equations
	)

	nButane := substance.NButane

	deviation := func(eos cubic.EOSType) float64 {
		cfg := nButane.CubicConfig(eos, zfactor.Args{T: T, P: measured, R: R})

		psat, err := cubic.SaturationPressure(cfg, T)
		if err != nil {
			t.Fatalf("SaturationPressure returned an unexpected error: %v", err)
		}

		return math.Abs(psat-measured) / measured
	}

	vdW := deviation(&cubic.VdW{})
	rk := deviation(&cubic.RK{})
	srk := deviation(&cubic.SRK{})
	pr := deviation(&cubic.PR{})

	if pr > acceptable {
		t.Errorf("Peng/Robinson is %.2f%% from the measured vapor pressure; want within %.0f%%", 100*pr, 100*acceptable)
	}

	if srk > acceptable {
		t.Errorf("Soave/Redlich/Kwong is %.2f%% from the measured vapor pressure; want within %.0f%%", 100*srk, 100*acceptable)
	}

	if !(pr < srk && srk < rk && rk < vdW) {
		t.Errorf(
			"expected accuracy to improve from vdW to PR; deviations were vdW %.2f%%, RK %.2f%%, SRK %.2f%%, PR %.2f%%",
			100*vdW, 100*rk, 100*srk, 100*pr,
		)
	}
}

// TestSaturationPressureAtCriticalPoint checks the boundary of the
// two-phase region.
//
// Above the critical temperature the phases are indistinguishable and no
// saturation pressure exists, so the critical pressure is returned
// instead of a failure.
func TestSaturationPressureAtCriticalPoint(t *testing.T) {
	nButane := substance.NButane
	cfg := nButane.CubicConfig(&cubic.PR{}, zfactor.Args{T: 350, P: 10, R: R})

	for _, temperature := range []float64{
		nButane.Critical.Tc,
		nButane.Critical.Tc + 1,
		nButane.Critical.Tc + 100,
	} {
		got, err := cubic.SaturationPressure(cfg, temperature)
		if err != nil {
			t.Fatalf("SaturationPressure returned an unexpected error: %v", err)
		}

		if got != nButane.Critical.Pc {
			t.Errorf("at %.2f K the saturation pressure is %.5f bar; want the critical pressure %.5f", temperature, got, nButane.Critical.Pc)
		}
	}
}

// TestSaturationPressureIncreasesWithTemperature checks that the
// predicted saturation curve rises monotonically, as the
// Clausius-Clapeyron relation requires.
func TestSaturationPressureIncreasesWithTemperature(t *testing.T) {
	nButane := substance.NButane
	cfg := nButane.CubicConfig(&cubic.PR{}, zfactor.Args{T: 350, P: 10, R: R})

	var previous float64

	for _, tr := range []float64{0.6, 0.7, 0.8, 0.9, 0.95} {
		temperature := tr * nButane.Critical.Tc

		got, err := cubic.SaturationPressure(cfg, temperature)
		if err != nil {
			t.Fatalf("SaturationPressure returned an unexpected error: %v", err)
		}

		if got <= previous {
			t.Errorf(
				"at Tr = %.2f the saturation pressure is %.6f bar; must exceed the previous value %.6f",
				tr, got, previous,
			)
		}

		if got >= nButane.Critical.Pc {
			t.Errorf("at Tr = %.2f the saturation pressure %.6f reaches the critical pressure", tr, got)
		}

		previous = got
	}
}

// TestSaturationPressureLowTemperatureLimit records where the iteration
// gives up, and checks that it does so cleanly.
//
// The pressure is refined by a multiplicative step damped to a factor of
// 1.2 either way, so descending to a saturation pressure many orders of
// magnitude below the critical pressure can exhaust the hundred
// iterations allowed. That happens for strongly non-ideal substances
// below roughly a third of their critical temperature, under the two
// temperature-dependent equations.
//
// What matters is not that convergence fails there — it is a region of
// little practical interest, and widening the damping would risk
// oscillation everywhere else — but that the failure is reported rather
// than dressed up as an answer.
func TestSaturationPressureLowTemperatureLimit(t *testing.T) {
	testCases := []struct {
		name string
		s    *substance.Substance
		eos  cubic.EOSType
		tr   float64
	}{
		{"ethanol/SRK", substance.Ethanol, &cubic.SRK{}, 0.35},
		{"ethanol/PR", substance.Ethanol, &cubic.PR{}, 0.35},
		{"n-octane/SRK", substance.NOctane, &cubic.SRK{}, 0.35},
		{"n-octane/PR", substance.NOctane, &cubic.PR{}, 0.35},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			temperature := tc.tr * tc.s.Critical.Tc
			cfg := tc.s.CubicConfig(tc.eos, zfactor.Args{T: temperature, P: tc.s.Critical.Pc, R: R})

			got, err := cubic.SaturationPressure(cfg, temperature)

			if err == nil {
				// Should the iteration be improved, the result must at
				// least be a physically sensible pressure.
				if got <= 0 || got >= tc.s.Critical.Pc {
					t.Errorf("converged to %g bar, which is outside (0, Pc)", got)
				}
				return
			}

			if got != 0 {
				t.Errorf("returned %g alongside the error; want 0", got)
			}
		})
	}
}

// TestLogFugacityIdealGasLimit checks that the fugacity coefficient
// approaches unity as the pressure falls.
//
// A gas becomes ideal in that limit, where fugacity and pressure
// coincide and their ratio's logarithm vanishes. The check covers van
// der Waals as well, whose equal sigma and epsilon send LogFugacity down
// a separate branch that the other equations never take.
func TestLogFugacityIdealGasLimit(t *testing.T) {
	const tol = 1e-4

	nButane := substance.NButane

	for name, eos := range map[string]cubic.EOSType{
		"vdW": &cubic.VdW{},
		"RK":  &cubic.RK{},
		"SRK": &cubic.SRK{},
		"PR":  &cubic.PR{},
	} {
		t.Run(name, func(t *testing.T) {
			const (
				temperature = 500.0
				pressure    = 1e-3
			)

			cfg := nButane.CubicConfig(eos, zfactor.Args{T: temperature, P: pressure, R: R})

			res, err := cubic.SolveForVolume(cfg)
			if err != nil {
				t.Fatalf("SolveForVolume returned an unexpected error: %v", err)
			}

			roots := res.Clean()

			RT := R * temperature
			a := res.A * pressure / (RT * RT)
			b := res.B * pressure / RT
			z := pressure * roots[len(roots)-1] / RT

			got := cubic.LogFugacity(cfg, z, a, b)

			if math.Abs(got) > tol {
				t.Errorf("ln(phi) = %.9f at %g bar; want approximately 0", got, pressure)
			}
		})
	}
}

// TestLogFugacityVanDerWaalsClosedForm checks the branch taken when the
// equation's sigma and epsilon coincide.
//
// The general expression divides by (epsilon - sigma), which vanishes for
// van der Waals, so LogFugacity substitutes the integral evaluated for
// that case. The result must match the closed form
//
//	ln(phi) = Z - 1 - ln(Z - B) - A/Z,
//
// which is the van der Waals fugacity coefficient.
func TestLogFugacityVanDerWaalsClosedForm(t *testing.T) {
	const tol = 1e-12

	nButane := substance.NButane

	for _, state := range []struct{ T, P float64 }{
		{350, 5}, {400, 20}, {500, 50},
	} {
		cfg := nButane.CubicConfig(&cubic.VdW{}, zfactor.Args{T: state.T, P: state.P, R: R})

		res, err := cubic.SolveForVolume(cfg)
		if err != nil {
			t.Fatalf("SolveForVolume returned an unexpected error: %v", err)
		}

		roots := res.Clean()

		RT := R * state.T
		a := res.A * state.P / (RT * RT)
		b := res.B * state.P / RT
		z := state.P * roots[len(roots)-1] / RT

		got := cubic.LogFugacity(cfg, z, a, b)
		want := z - 1 - math.Log(z-b) - a/z

		if math.Abs(got-want) > tol {
			t.Errorf(
				"at %g K and %g bar: ln(phi) = %.12f; the van der Waals closed form gives %.12f",
				state.T, state.P, got, want,
			)
		}
	}
}

// TestLogFugacityBelowUnityWhereAttractionDominates checks the sign of
// the fugacity coefficient for a gas.
//
// Attractive intermolecular forces lower the fugacity below the
// pressure, so the coefficient is less than one and its logarithm
// negative wherever those forces dominate. A sign error in the
// integrated term would show here while leaving the magnitude plausible.
func TestLogFugacityBelowUnityWhereAttractionDominates(t *testing.T) {
	nButane := substance.NButane

	for name, eos := range map[string]cubic.EOSType{
		"vdW": &cubic.VdW{},
		"RK":  &cubic.RK{},
		"SRK": &cubic.SRK{},
		"PR":  &cubic.PR{},
	} {
		t.Run(name, func(t *testing.T) {
			const (
				temperature = 450.0
				pressure    = 10.0
			)

			cfg := nButane.CubicConfig(eos, zfactor.Args{T: temperature, P: pressure, R: R})

			res, err := cubic.SolveForVolume(cfg)
			if err != nil {
				t.Fatalf("SolveForVolume returned an unexpected error: %v", err)
			}

			roots := res.Clean()

			RT := R * temperature
			a := res.A * pressure / (RT * RT)
			b := res.B * pressure / RT
			z := pressure * roots[len(roots)-1] / RT

			if got := cubic.LogFugacity(cfg, z, a, b); got >= 0 {
				t.Errorf("ln(phi) = %.9f; want a negative value where attraction dominates", got)
			}
		})
	}
}

// TestConfigConstructors checks that each convenience constructor
// records the conditions and critical properties it is given, and
// selects the matching equation of state.
//
// The two equations without a temperature-dependent alpha take no
// acentric factor, so theirs is left at zero; nothing reads it.
func TestConfigConstructors(t *testing.T) {
	const (
		T  = 350.0
		P  = 9.4573
		Tc = 425.1
		Pc = 37.96
		w  = 0.2
	)

	testCases := []struct {
		name         string
		cfg          *cubic.EOSCfg
		wantType     any
		wantAcentric float64
	}{
		{"vdW", cubic.NewvdWCfg(T, P, Tc, Pc, R), &cubic.VdW{}, 0},
		{"RK", cubic.NewRKCfg(T, P, Tc, Pc, R), &cubic.RK{}, 0},
		{"SRK", cubic.NewSRKCfg(T, P, Tc, Pc, w, R), &cubic.SRK{}, w},
		{"PR", cubic.NewPRCfg(T, P, Tc, Pc, w, R), &cubic.PR{}, w},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if got, want := fmt.Sprintf("%T", tc.cfg.Type), fmt.Sprintf("%T", tc.wantType); got != want {
				t.Errorf("Type = %s; want %s", got, want)
			}

			if tc.cfg.T != T || tc.cfg.P != P {
				t.Errorf("conditions = (%g, %g); want (%g, %g)", tc.cfg.T, tc.cfg.P, T, P)
			}

			if tc.cfg.Tc != Tc || tc.cfg.Pc != Pc {
				t.Errorf("critical properties = (%g, %g); want (%g, %g)", tc.cfg.Tc, tc.cfg.Pc, Tc, Pc)
			}

			if tc.cfg.R != R {
				t.Errorf("R = %g; want %g", tc.cfg.R, R)
			}

			if tc.cfg.Acentric != tc.wantAcentric {
				t.Errorf("acentric factor = %g; want %g", tc.cfg.Acentric, tc.wantAcentric)
			}

			// A configuration from a constructor must be usable.
			if _, err := cubic.SolveForVolume(tc.cfg); err != nil {
				t.Errorf("SolveForVolume returned an unexpected error: %v", err)
			}
		})
	}
}

// TestResultStringers checks that the result types render their values,
// so that a printed result is informative rather than an address.
func TestResultStringers(t *testing.T) {
	cfg := substance.NButane.CubicConfig(&cubic.PR{}, zfactor.Args{T: 350, P: 9.4573, R: R})

	volume, err := cubic.SolveForVolume(cfg)
	if err != nil {
		t.Fatalf("SolveForVolume returned an unexpected error: %v", err)
	}

	if got := volume.String(); !containsAll(got, "VolumeResult", "A:", "B:", "Volumes:") {
		t.Errorf("VolumeResult.String() = %q; want the parameters and roots", got)
	}

	pressure, err := cubic.Pressure(cfg, volume.Clean()[0])
	if err != nil {
		t.Fatalf("Pressure returned an unexpected error: %v", err)
	}

	if got := pressure.String(); !containsAll(got, "PressureResult", "A:", "B:", "P:") {
		t.Errorf("PressureResult.String() = %q; want the parameters and pressure", got)
	}
}

// TestPressureInvalidInput checks the guards on the pressure
// calculation, which shares its validation with SolveForVolume but is
// reached separately.
func TestPressureInvalidInput(t *testing.T) {
	testCases := []struct {
		name string
		cfg  *cubic.EOSCfg
	}{
		{"non-positive temperature", &cubic.EOSCfg{Type: &cubic.PR{}, T: 0, P: 10, Tc: 425.1, Pc: 37.96, R: R}},
		{"invalid critical temperature", &cubic.EOSCfg{Type: &cubic.PR{}, T: 350, P: 10, Tc: 0, Pc: 37.96, R: R}},
		{"invalid critical pressure", &cubic.EOSCfg{Type: &cubic.PR{}, T: 350, P: 10, Tc: 425.1, Pc: 0, R: R}},
		{"invalid gas constant", &cubic.EOSCfg{Type: &cubic.PR{}, T: 350, P: 10, Tc: 425.1, Pc: 37.96, R: 0}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := cubic.Pressure(tc.cfg, 200); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// containsAll reports whether s contains every one of the substrings.
func containsAll(s string, substrings ...string) bool {
	for _, sub := range substrings {
		if !strings.Contains(s, sub) {
			return false
		}
	}
	return true
}
