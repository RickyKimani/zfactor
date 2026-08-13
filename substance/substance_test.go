package substance_test

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
	"github.com/rickykimani/zfactor/substance"
)

// R is the gas constant in bar·cm³/(mol·K), consistent with the units
// used by the substance table: critical pressures in bar, critical
// volumes in cm³/mol and temperatures in kelvin.
const R = 10 * zfactor.RSI

// TestAllRegistry checks the generated registry is usable: populated,
// free of nil entries, and free of duplicate names.
//
// Duplicate names matter because the generator derives Go identifiers
// from them, so a collision would either fail to compile or silently
// shadow an entry.
func TestAllRegistry(t *testing.T) {
	all := substance.All

	if len(all) == 0 {
		t.Fatal("substance.All is empty")
	}

	seen := make(map[string]int, len(all))

	for i, s := range all {
		if s == nil {
			t.Fatalf("substance.All[%d] is nil", i)
		}

		if s.Name == "" {
			t.Errorf("substance.All[%d] has an empty name", i)
			continue
		}

		if first, ok := seen[s.Name]; ok {
			t.Errorf("duplicate substance name %q at indices %d and %d", s.Name, first, i)
			continue
		}

		seen[s.Name] = i
	}
}

// TestCriticalPropertiesPositive checks that every tabulated property
// is physically meaningful.
//
// A zero or negative value here would propagate silently: the
// equation-of-state and correlation code divides by Tc and Pc when
// forming reduced conditions.
func TestCriticalPropertiesPositive(t *testing.T) {
	for _, s := range substance.All {
		t.Run(s.Name, func(t *testing.T) {
			properties := []struct {
				name  string
				value float64
			}{
				{"MW", s.MW},
				{"Tc", s.Critical.Tc},
				{"Pc", s.Critical.Pc},
				{"Vc", s.Critical.Vc},
				{"Zc", s.Critical.Zc},
			}

			for _, p := range properties {
				if p.value <= 0 {
					t.Errorf("%s = %g; want a positive value", p.name, p.value)
				}
			}
		})
	}
}

// zcExceptions lists substances whose tabulated critical compressibility
// factor is inconsistent with their own critical constants in the source
// data.
//
// Appendix B gives Zc = 0.265 for chlorine, which implies a critical
// volume of 119.2 cm³/mol rather than the tabulated 124.0. The two
// values evidently come from different determinations; the critical
// volume is the least precisely known of the critical properties. The
// entry is reproduced faithfully from the source and excluded here
// rather than treated as a transcription error.
var zcExceptions = map[string]string{
	"Chlorine": "Zc and Vc are mutually inconsistent in Appendix B",
}

// TestCompressibilityFactorConsistency checks each tabulated Zc against
// its definition,
//
//	Zc = Pc Vc / (R Tc),
//
// which is an identity rather than a correlation.
//
// The purpose is to guard the data pipeline rather than the physics.
// Should the generator ever mismap a column, apply the wrong unit scale
// or slip by one row, every entry would break by a wide margin. The
// tolerance is therefore loose enough to accommodate the rounding and
// mixed sourcing present in the published table — all but one entry
// agree to within 1% — while remaining far tighter than any plausible
// pipeline defect.
func TestCompressibilityFactorConsistency(t *testing.T) {
	const tol = 0.02

	for _, s := range substance.All {
		t.Run(s.Name, func(t *testing.T) {
			if reason, skip := zcExceptions[s.Name]; skip {
				t.Skipf("known source inconsistency: %s", reason)
			}

			c := s.Critical
			want := c.Pc * c.Vc / (R * c.Tc)

			if got := math.Abs(want-c.Zc) / c.Zc; got > tol {
				t.Errorf(
					"Zc = %.4f but Pc*Vc/(R*Tc) = %.4f (%.2f%% apart)",
					c.Zc, want, 100*got,
				)
			}
		})
	}
}

// TestAcentricFactorRange checks that acentric factors fall within the
// range spanned by real fluids.
//
// The scale is defined so that the simple fluids (argon, krypton,
// xenon) sit at zero. Quantum fluids fall below it — helium is the
// lowest entry in the table — and strongly polar or associating fluids
// rise above it, but values outside these bounds indicate corrupt data.
func TestAcentricFactorRange(t *testing.T) {
	const (
		low  = -0.5
		high = 1.0
	)

	for _, s := range substance.All {
		if s.Acentric < low || s.Acentric > high {
			t.Errorf(
				"%s: acentric factor = %.4f; want between %.1f and %.1f",
				s.Name, s.Acentric, low, high,
			)
		}
	}
}

// TestNormalBoilingPointBelowCritical checks that tabulated normal
// boiling points lie below the critical temperature, as they must: a
// liquid cannot boil above the critical point.
//
// Entries with no tabulated boiling point carry zero and are skipped.
func TestNormalBoilingPointBelowCritical(t *testing.T) {
	for _, s := range substance.All {
		if s.Tn == 0 {
			continue
		}

		if s.Tn >= s.Critical.Tc {
			t.Errorf(
				"%s: normal boiling point %.2f K is not below the critical temperature %.2f K",
				s.Name, s.Tn, s.Critical.Tc,
			)
		}
	}
}

// TestVaporPressureAtNormalBoilingPoint exercises the Lee-Kesler vapor
// pressure correlation against its own definition.
//
// The normal boiling point is the temperature at which the saturation
// pressure equals one atmosphere, and the correlation estimates the
// acentric factor from that same point. Evaluating it at Tn must
// therefore return one atmosphere exactly, for every substance and
// independently of how well the correlation performs elsewhere.
func TestVaporPressureAtNormalBoilingPoint(t *testing.T) {
	const tol = 1e-9

	for _, s := range substance.All {
		if s.Tn == 0 {
			continue
		}

		t.Run(s.Name, func(t *testing.T) {
			got, err := s.LeeKeslerVaporPressure(s.Tn)
			if err != nil {
				t.Fatalf("LeeKeslerVaporPressure returned an unexpected error: %v", err)
			}

			if rel := math.Abs(got-zfactor.AtmBar) / zfactor.AtmBar; rel > tol {
				t.Errorf(
					"vapor pressure at the normal boiling point = %.8f bar; want %.8f bar (1 atm)",
					got, zfactor.AtmBar,
				)
			}
		})
	}
}

// TestVsatAtCriticalTemperature exercises the Rackett equation against
// its own limiting case.
//
// Rackett gives Vsat = Vc·Zc^((1-Tr)^(2/7)), so at the critical
// temperature the exponent vanishes and the saturated liquid volume
// must equal the critical volume exactly.
func TestVsatAtCriticalTemperature(t *testing.T) {
	const tol = 1e-12

	for _, s := range substance.All {
		t.Run(s.Name, func(t *testing.T) {
			got, err := s.Vsat(s.Critical.Tc)
			if err != nil {
				t.Fatalf("Vsat returned an unexpected error: %v", err)
			}

			if rel := math.Abs(got-s.Critical.Vc) / s.Critical.Vc; rel > tol {
				t.Errorf(
					"Vsat(Tc) = %.6f cm³/mol; want the critical volume %.6f",
					got, s.Critical.Vc,
				)
			}
		})
	}
}

// TestVsatBelowCriticalVolume checks the physical ordering of the
// Rackett equation: a saturated liquid is denser than the critical
// fluid, so its molar volume must be smaller, and it expands as the
// temperature rises toward the critical point.
func TestVsatBelowCriticalVolume(t *testing.T) {
	for _, s := range substance.All {
		t.Run(s.Name, func(t *testing.T) {
			var previous float64

			for _, tr := range []float64{0.5, 0.6, 0.7, 0.8, 0.9} {
				got, err := s.Vsat(tr * s.Critical.Tc)
				if err != nil {
					t.Fatalf("Vsat returned an unexpected error: %v", err)
				}

				if got >= s.Critical.Vc {
					t.Errorf(
						"Vsat at Tr = %.1f is %.4f cm³/mol; must be below the critical volume %.4f",
						tr, got, s.Critical.Vc,
					)
				}

				if got <= previous {
					t.Errorf(
						"Vsat at Tr = %.1f is %.4f cm³/mol; must exceed the previous value %.4f",
						tr, got, previous,
					)
				}

				previous = got
			}
		})
	}
}

// TestLeeKeslerAcentricEstimate checks the acentric factor estimated
// from the normal boiling point against the tabulated value.
//
// The correlation is fitted to normal fluids and is not expected to
// reproduce strongly associating ones — sulfuric acid, the worst entry
// in the table, is off by roughly 0.48 — so only well-behaved
// substances are checked here.
func TestLeeKeslerAcentricEstimate(t *testing.T) {
	const tol = 0.05

	testCases := []*substance.Substance{
		substance.Methane,
		substance.Ethane,
		substance.Propane,
		substance.NButane,
		substance.CarbonDioxide,
	}

	for _, s := range testCases {
		t.Run(s.Name, func(t *testing.T) {
			if s.Tn == 0 {
				t.Skip("no tabulated normal boiling point")
			}

			got, err := s.LeeKeslerAcentric()
			if err != nil {
				t.Fatalf("LeeKeslerAcentric returned an unexpected error: %v", err)
			}

			if math.Abs(got-s.Acentric) > tol {
				t.Errorf(
					"estimated acentric factor = %.4f; tabulated value is %.4f",
					got, s.Acentric,
				)
			}
		})
	}
}

// TestLeeKeslerExample3_10 reproduces part (b) of Example 3.10 of Smith,
// Van Ness & Abbott: the molar volume of n-butane at 510 K and 25 bar
// from the generalized compressibility-factor correlation.
//
// Unlike the test in the lee-kesler package, which supplies the reduced
// conditions directly, this exercises the whole path: the substance's
// own critical constants reduce the state, the correlation is evaluated,
// and the tabulated acentric factor weights the departure term. The
// example rounds Tr and Pr to three decimals before reading the tables,
// so the tolerance accommodates that rather than the interpolation.
func TestLeeKeslerExample3_10(t *testing.T) {
	const (
		T = 510.0 // K
		P = 25.0  // bar

		wantZ = 0.873
		wantV = 1480.7 // cm³/mol

		zTol   = 2e-3
		relTol = 2e-3
	)

	got, err := substance.NButane.LeeKesler(
		zfactor.Args{T: T, P: P},
		leekesler.CompressibilityFactor,
	)
	if err != nil {
		t.Fatalf("LeeKesler returned an unexpected error: %v", err)
	}

	if math.Abs(got-wantZ) > zTol {
		t.Errorf("compressibility factor = %.6f; want %.3f", got, wantZ)
	}

	v := got * R * T / P

	if rel := math.Abs(v-wantV) / wantV; rel > relTol {
		t.Errorf("molar volume = %.2f cm³/mol; want %.1f (%.3f%% apart)", v, wantV, 100*rel)
	}
}

// TestMissingNormalBoilingPoint checks that the correlations requiring a
// normal boiling point refuse to run without one, rather than treating
// the absent value as a temperature of zero.
func TestMissingNormalBoilingPoint(t *testing.T) {
	noBoilingPoint := &substance.Substance{
		Name:     "test substance",
		MW:       50,
		Acentric: 0.1,
		Critical: substance.CriticalProps{Tc: 400, Pc: 40, Vc: 200, Zc: 0.27},
	}

	if _, err := noBoilingPoint.LeeKeslerAcentric(); err == nil {
		t.Error("LeeKeslerAcentric: expected an error; got nil")
	}

	if _, err := noBoilingPoint.LeeKeslerVaporPressure(350); err == nil {
		t.Error("LeeKeslerVaporPressure: expected an error; got nil")
	}
}

// TestInvalidConditions checks that non-physical states are rejected.
func TestInvalidConditions(t *testing.T) {
	s := substance.NButane

	t.Run("Vsat with non-positive temperature", func(t *testing.T) {
		if _, err := s.Vsat(0); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("ReducedDensity with non-positive temperature", func(t *testing.T) {
		if _, err := s.ReducedDensity(zfactor.Args{T: 0, P: 10}); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("ReducedDensity with negative pressure", func(t *testing.T) {
		if _, err := s.ReducedDensity(zfactor.Args{T: 350, P: -1}); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("AbbottResidualEnthalpy with non-positive temperature", func(t *testing.T) {
		if _, err := s.AbbottResidualEnthalpy(zfactor.Args{T: 0, P: 10}); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("AbbottResidualEntropy with non-positive pressure", func(t *testing.T) {
		if _, err := s.AbbottResidualEntropy(zfactor.Args{T: 350, P: 0}); err == nil {
			t.Error("expected an error; got nil")
		}
	})
}

// TestCubicConfig checks that the configuration carries the substance's
// properties and the supplied conditions through to the solver, for
// every equation of state.
func TestCubicConfig(t *testing.T) {
	s := substance.NButane
	args := zfactor.Args{T: 350, P: 9.4573, R: R}

	for _, eos := range []cubic.EOSType{
		&cubic.VdW{}, &cubic.RK{}, &cubic.SRK{}, &cubic.PR{},
	} {
		cfg := s.CubicConfig(eos, args)

		if cfg.Type != eos {
			t.Errorf("Type = %T; want %T", cfg.Type, eos)
		}

		if cfg.T != args.T || cfg.P != args.P || cfg.R != args.R {
			t.Errorf(
				"conditions = (T %g, P %g, R %g); want (%g, %g, %g)",
				cfg.T, cfg.P, cfg.R, args.T, args.P, args.R,
			)
		}

		if cfg.Tc != s.Critical.Tc || cfg.Pc != s.Critical.Pc {
			t.Errorf(
				"critical properties = (Tc %g, Pc %g); want (%g, %g)",
				cfg.Tc, cfg.Pc, s.Critical.Tc, s.Critical.Pc,
			)
		}

		if cfg.Acentric != s.Acentric {
			t.Errorf("acentric factor = %g; want %g", cfg.Acentric, s.Acentric)
		}
	}
}

// TestNewLinearMixtureSingleComponent checks the limiting case of Kay's
// rule: a mixture consisting of one substance at unit mole fraction
// must reproduce that substance's properties exactly.
func TestNewLinearMixtureSingleComponent(t *testing.T) {
	want := substance.Propane

	got, err := substance.NewLinearMixture("pure propane", []substance.Component{
		{Substance: want, Fraction: 1.0},
	})
	if err != nil {
		t.Fatalf("NewLinearMixture returned an unexpected error: %v", err)
	}

	if got.MW != want.MW {
		t.Errorf("MW = %g; want %g", got.MW, want.MW)
	}

	if got.Acentric != want.Acentric {
		t.Errorf("acentric factor = %g; want %g", got.Acentric, want.Acentric)
	}

	if got.Critical != want.Critical {
		t.Errorf("critical properties = %+v; want %+v", got.Critical, want.Critical)
	}
}

// TestNewLinearMixtureAverages checks the pseudo-critical properties of
// an equimolar mixture against the linear averages Kay's rule defines.
func TestNewLinearMixtureAverages(t *testing.T) {
	a, b := substance.CarbonDioxide, substance.Propane

	got, err := substance.NewLinearMixture("equimolar CO2/propane", []substance.Component{
		{Substance: a, Fraction: 0.5},
		{Substance: b, Fraction: 0.5},
	})
	if err != nil {
		t.Fatalf("NewLinearMixture returned an unexpected error: %v", err)
	}

	const tol = 1e-9

	testCases := []struct {
		name      string
		got, want float64
	}{
		{"MW", got.MW, (a.MW + b.MW) / 2},
		{"acentric factor", got.Acentric, (a.Acentric + b.Acentric) / 2},
		{"Tc", got.Critical.Tc, (a.Critical.Tc + b.Critical.Tc) / 2},
		{"Pc", got.Critical.Pc, (a.Critical.Pc + b.Critical.Pc) / 2},
		{"Vc", got.Critical.Vc, (a.Critical.Vc + b.Critical.Vc) / 2},
		{"Zc", got.Critical.Zc, (a.Critical.Zc + b.Critical.Zc) / 2},
	}

	for _, tc := range testCases {
		if math.Abs(tc.got-tc.want) > tol {
			t.Errorf("%s = %g; want %g", tc.name, tc.got, tc.want)
		}
	}
}

// TestNewLinearMixtureInvalid checks that malformed mixtures are
// rejected.
func TestNewLinearMixtureInvalid(t *testing.T) {
	testCases := []struct {
		name       string
		components []substance.Component
	}{
		{
			name:       "no components",
			components: nil,
		},
		{
			name:       "nil substance",
			components: []substance.Component{{Substance: nil, Fraction: 1.0}},
		},
		{
			name: "negative mole fraction",
			components: []substance.Component{
				{Substance: substance.Methane, Fraction: -0.5},
				{Substance: substance.Ethane, Fraction: 1.5},
			},
		},
		{
			name: "mole fractions do not sum to one",
			components: []substance.Component{
				{Substance: substance.Methane, Fraction: 0.3},
				{Substance: substance.Ethane, Fraction: 0.3},
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := substance.NewLinearMixture("invalid", tc.components); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}
