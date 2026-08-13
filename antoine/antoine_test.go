package antoine_test

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/antoine"
)

// atmKPa is one standard atmosphere in kPa, the pressure that defines
// the normal boiling point.
const atmKPa = zfactor.AtmKPa

// models returns a representative spread of the generated table:
// substances covering the full range of volatility, plus every entry
// whose formula ends in a subscript, since those are the ones a parsing
// defect reaches first.
func models() map[string]*antoine.Antoine {
	return map[string]*antoine.Antoine{
		"Water":           antoine.Water,
		"Ethanol":         antoine.Ethanol,
		"Benzene":         antoine.Benzene,
		"Acetone":         antoine.Acetone,
		"Methanol":        antoine.Methanol,
		"Toluene":         antoine.Toluene,
		"AceticAcid":      antoine.AceticAcid,
		"Dichloromethane": antoine.Dichloromethane,
		"MethylAcetate":   antoine.MethylAcetate,
		"Nitromethane":    antoine.Nitromethane,
		"EthyleneGlycol":  antoine.EthyleneGlycol,
		"One4Dioxane":     antoine.One4Dioxane,
	}
}

// TestRoundTrip checks that Temperature inverts Pressure exactly.
//
// The two are closed-form inverses of the same equation, so composing
// them must return the original temperature to within rounding. This
// holds for any substance and any temperature, independently of how well
// the correlation reproduces experiment.
func TestRoundTrip(t *testing.T) {
	const tol = 1e-9

	for name, model := range models() {
		t.Run(name, func(t *testing.T) {
			span := model.Range.High - model.Range.Low

			for _, fraction := range []float64{0, 0.1, 0.25, 0.5, 0.75, 0.9, 1} {
				want := model.Range.Low + fraction*span

				p, err := model.Pressure(want)
				if err != nil {
					t.Fatalf("Pressure(%g) returned an unexpected error: %v", want, err)
				}

				got, err := model.Temperature(p)
				if err != nil {
					t.Fatalf("Temperature(%g) returned an unexpected error: %v", p, err)
				}

				if math.Abs(got-want) > tol {
					t.Errorf(
						"Temperature(Pressure(%g)) = %.12f; want %.12f",
						want, got, want,
					)
				}
			}
		})
	}
}

// TestNormalBoilingPoint checks each correlation against its own
// tabulated normal boiling point.
//
// The normal boiling point is by definition the temperature at which the
// saturation pressure equals one atmosphere, so evaluating the
// correlation there must return 101.325 kPa. The coefficients and the
// boiling point are separate columns of the source table, so nothing but
// this check ties them together.
//
// It is not a formality. Six entries had a trailing formula subscript
// absorbed into their A coefficient — "C2H4O2" followed by "15.0717"
// extracts from the PDF as the single token "C2H4O215.0717" — which left
// A too large by 200 and the predicted pressure wrong by a factor of
// 10^87. This check identifies every one of them.
func TestNormalBoilingPoint(t *testing.T) {
	const relTol = 5e-3

	for name, model := range models() {
		t.Run(name, func(t *testing.T) {
			if model.Tn == 0 {
				t.Skip("no tabulated normal boiling point")
			}

			got, err := model.Pressure(model.Tn)

			// The boiling point may sit outside the fitted range, which
			// is a caveat on accuracy rather than a failure.
			var rangeErr *antoine.RangeError
			if err != nil && !errors.As(err, &rangeErr) {
				t.Fatalf("Pressure returned an unexpected error: %v", err)
			}

			if rel := math.Abs(got-atmKPa) / atmKPa; rel > relTol {
				t.Errorf(
					"pressure at the normal boiling point of %g °C = %.4f kPa; want %.4f kPa (%.3f%% apart)",
					model.Tn, got, atmKPa, 100*rel,
				)
			}
		})
	}
}

// TestCoefficientsAreSane checks that the constants fall in the range
// every substance in the table occupies.
//
// In this form of the equation, ln(P/kPa) = A - B/(t + C), A is the
// logarithm of the pressure the substance would exert at infinite
// temperature and sits between roughly 12 and 20; B and C are both
// positive. A value far outside those bounds is not a different
// substance but corrupt data.
func TestCoefficientsAreSane(t *testing.T) {
	for name, model := range models() {
		t.Run(name, func(t *testing.T) {
			if model.A < 10 || model.A > 25 {
				t.Errorf("A = %g; want a value between 10 and 25", model.A)
			}

			if model.B <= 0 {
				t.Errorf("B = %g; want a positive value", model.B)
			}

			if model.C <= 0 {
				t.Errorf("C = %g; want a positive value", model.C)
			}

			if model.Range.High <= model.Range.Low {
				t.Errorf(
					"range [%g, %g] is empty",
					model.Range.Low, model.Range.High,
				)
			}

			// The pole of the correlation sits at t = -C, which must lie
			// below the fitted range for the equation to be usable
			// across it.
			if -model.C >= model.Range.Low {
				t.Errorf(
					"the pole at t = %g lies inside the fitted range [%g, %g]",
					-model.C, model.Range.Low, model.Range.High,
				)
			}
		})
	}
}

// TestFormulasKeepTheirSubscripts pins the formulas of the entries whose
// trailing subscript was previously absorbed into the A coefficient.
//
// The formula is only a label, but it is the visible half of the same
// defect TestNormalBoilingPoint detects numerically, and it is the
// cheaper of the two to read when something has gone wrong.
func TestFormulasKeepTheirSubscripts(t *testing.T) {
	testCases := []struct {
		model *antoine.Antoine
		want  string
	}{
		{antoine.AceticAcid, "C2H4O2"},
		{antoine.Dichloromethane, "CH2Cl2"},
		{antoine.MethylAcetate, "C3H6O2"},
		{antoine.Nitromethane, "CH3NO2"},
		{antoine.EthyleneGlycol, "C2H6O2"},
		{antoine.One4Dioxane, "C4H8O2"},
		{antoine.Water, "H2O"},
		{antoine.Benzene, "C6H6"},
	}

	for _, tc := range testCases {
		t.Run(tc.want, func(t *testing.T) {
			if got := tc.model.Formula; got != tc.want {
				t.Errorf("formula = %q; want %q", got, tc.want)
			}
		})
	}
}

// TestPressureIncreasesWithTemperature checks that each saturation curve
// rises monotonically, as the Clausius-Clapeyron relation requires of
// any substance.
func TestPressureIncreasesWithTemperature(t *testing.T) {
	for name, model := range models() {
		t.Run(name, func(t *testing.T) {
			span := model.Range.High - model.Range.Low

			var previous float64

			for i := range 11 {
				temperature := model.Range.Low + span*float64(i)/10

				got, err := model.Pressure(temperature)
				if err != nil {
					t.Fatalf("Pressure returned an unexpected error: %v", err)
				}

				if got <= previous {
					t.Errorf(
						"at %g °C the pressure is %.6f kPa; must exceed the previous value %.6f",
						temperature, got, previous,
					)
				}

				previous = got
			}
		})
	}
}

// TestRangeErrorAccompaniesResult checks the contract for temperatures
// outside the fitted interval: the value is returned together with a
// *RangeError rather than discarded.
//
// The correlation stays defined outside its range and callers may
// legitimately accept a small extrapolation — the VLE solvers do exactly
// that, treating a RangeError as a caveat and any other error as fatal.
func TestRangeErrorAccompaniesResult(t *testing.T) {
	model := antoine.Water

	for _, temperature := range []float64{model.Range.Low - 5, model.Range.High + 5} {
		got, err := model.Pressure(temperature)

		if err == nil {
			t.Fatalf("Pressure(%g): expected a range error; got nil", temperature)
		}

		var rangeErr *antoine.RangeError
		if !errors.As(err, &rangeErr) {
			t.Fatalf("Pressure(%g): error is not a *RangeError: %v", temperature, err)
		}

		if rangeErr.T != temperature {
			t.Errorf("RangeError reports T = %g; want %g", rangeErr.T, temperature)
		}

		if rangeErr.Low != model.Range.Low || rangeErr.High != model.Range.High {
			t.Errorf(
				"RangeError reports [%g, %g]; want [%g, %g]",
				rangeErr.Low, rangeErr.High, model.Range.Low, model.Range.High,
			)
		}

		if got <= 0 {
			t.Errorf("Pressure(%g) = %g; want the extrapolated value alongside the error", temperature, got)
		}
	}
}

// TestValidateTempRange checks the bounds are inclusive.
func TestValidateTempRange(t *testing.T) {
	model := antoine.Water

	testCases := []struct {
		name string
		t    float64
		want bool
	}{
		{"below the range", model.Range.Low - 1, false},
		{"at the lower bound", model.Range.Low, true},
		{"inside the range", (model.Range.Low + model.Range.High) / 2, true},
		{"at the upper bound", model.Range.High, true},
		{"above the range", model.Range.High + 1, false},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if got := model.ValidateTempRange(tc.t); got != tc.want {
				t.Errorf("ValidateTempRange(%g) = %v; want %v", tc.t, got, tc.want)
			}
		})
	}
}

// TestTemperatureInvalidPressure checks that pressures with no
// corresponding saturation temperature are refused.
//
// Below zero the logarithm is undefined. At and beyond exp(A) the
// denominator of the inverted equation vanishes and then changes sign,
// which would otherwise yield an infinity or a temperature below the
// pole — a plausible-looking number describing nothing.
func TestTemperatureInvalidPressure(t *testing.T) {
	model := antoine.Water

	limit := math.Exp(model.A)

	testCases := []struct {
		name string
		p    float64
	}{
		{"zero", 0},
		{"negative", -10},
		{"at the limit", limit},
		{"beyond the limit", limit * 10},
		{"absurdly large", 1e300},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := model.Temperature(tc.p)

			if err == nil {
				t.Fatalf("Temperature(%g) = %g; expected an error", tc.p, got)
			}

			if got != 0 {
				t.Errorf("Temperature(%g) = %g; want 0 alongside the error", tc.p, got)
			}
		})
	}
}

// TestKnownVaporPressures checks a few widely quoted saturation
// pressures, so the table is anchored to values independent of the
// correlation's own internal consistency.
func TestKnownVaporPressures(t *testing.T) {
	testCases := []struct {
		name        string
		model       *antoine.Antoine
		temperature float64
		want        float64
		relTol      float64
	}{
		// Water boils at 100 °C under one atmosphere.
		{"water at 100 C", antoine.Water, 100, atmKPa, 5e-3},
		// Water at ambient temperature, roughly 3.17 kPa.
		{"water at 25 C", antoine.Water, 25, 3.17, 2e-2},
		// Ethanol at ambient temperature, roughly 7.9 kPa.
		{"ethanol at 25 C", antoine.Ethanol, 25, 7.9, 3e-2},
		// Benzene at ambient temperature, roughly 12.7 kPa.
		{"benzene at 25 C", antoine.Benzene, 25, 12.7, 3e-2},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := tc.model.Pressure(tc.temperature)
			if err != nil {
				t.Fatalf("Pressure returned an unexpected error: %v", err)
			}

			if rel := math.Abs(got-tc.want) / tc.want; rel > tc.relTol {
				t.Errorf(
					"pressure = %.4f kPa; want approximately %.4f (%.2f%% apart)",
					got, tc.want, 100*rel,
				)
			}
		})
	}
}

// TestModelInterface checks that *Antoine satisfies the Model interface
// the VLE packages accept.
func TestModelInterface(t *testing.T) {
	var model antoine.Model = antoine.Water

	if _, err := model.Pressure(50); err != nil {
		t.Errorf("Pressure returned an unexpected error: %v", err)
	}

	if _, err := model.LnPSat(50); err != nil {
		t.Errorf("LnPSat returned an unexpected error: %v", err)
	}

	if _, err := model.Temperature(atmKPa); err != nil {
		t.Errorf("Temperature returned an unexpected error: %v", err)
	}

	if !model.ValidateTempRange(50) {
		t.Error("ValidateTempRange(50) = false; want true for water")
	}
}

// TestLnPSatMatchesPressure checks the two accessors agree, since
// Pressure is defined as the exponential of LnPSat.
func TestLnPSatMatchesPressure(t *testing.T) {
	const tol = 1e-12

	for name, model := range models() {
		t.Run(name, func(t *testing.T) {
			temperature := (model.Range.Low + model.Range.High) / 2

			lnP, err := model.LnPSat(temperature)
			if err != nil {
				t.Fatalf("LnPSat returned an unexpected error: %v", err)
			}

			p, err := model.Pressure(temperature)
			if err != nil {
				t.Fatalf("Pressure returned an unexpected error: %v", err)
			}

			if got := math.Exp(lnP); math.Abs(got-p) > tol*math.Max(p, 1) {
				t.Errorf("exp(LnPSat) = %.12f but Pressure = %.12f", got, p)
			}
		})
	}
}
