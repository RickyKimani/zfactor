package liquids_test

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor/liquids"
)

// TestVsatAtCriticalTemperature checks the Rackett equation against its
// own limiting case.
//
// The correlation is Vsat = Vc·Zc^((1-Tr)^(2/7)), so at the critical
// temperature the exponent vanishes, Zc is raised to the zeroth power
// and the saturated liquid volume equals the critical volume. This holds
// for any substance and is exact.
func TestVsatAtCriticalTemperature(t *testing.T) {
	const tol = 1e-12

	testCases := []struct {
		name   string
		vc, zc float64
	}{
		{"n-butane", 255.0, 0.274},
		{"carbon dioxide", 94.0, 0.274},
		{"water", 55.9, 0.229},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := liquids.Vsat(tc.vc, tc.zc, 1)
			if err != nil {
				t.Fatalf("Vsat returned an unexpected error: %v", err)
			}

			if math.Abs(got-tc.vc) > tol {
				t.Errorf("Vsat at Tr = 1 is %.9f; want the critical volume %.9f", got, tc.vc)
			}
		})
	}
}

// TestVsatBelowCriticalVolume checks that a saturated liquid is denser
// than the critical fluid and expands as it warms toward the critical
// point.
//
// Zc is below one for every real substance, so raising it to a positive
// power gives a factor below one, and that power falls as Tr rises.
func TestVsatIncreasesTowardTheCriticalVolume(t *testing.T) {
	const (
		vc = 255.0
		zc = 0.274
	)

	var previous float64

	for _, tr := range []float64{0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95} {
		got, err := liquids.Vsat(vc, zc, tr)
		if err != nil {
			t.Fatalf("Vsat returned an unexpected error: %v", err)
		}

		if got >= vc {
			t.Errorf("at Tr = %.2f the volume %.4f reaches the critical volume %.4f", tr, got, vc)
		}

		if got <= previous {
			t.Errorf("at Tr = %.2f the volume %.4f does not exceed the previous %.4f", tr, got, previous)
		}

		previous = got
	}
}

// TestVsatAboveCriticalTemperature checks the squaring that keeps the
// correlation defined past the critical point.
//
// The exponent (1-Tr)^(2/7) has no real value for Tr above one, since a
// negative base is raised to a fractional power. The implementation
// squares the base first, which is equivalent below the critical
// temperature and finite above it. The result no longer describes a
// saturated liquid there, but it must at least be a positive number
// rather than a NaN, and symmetric about Tr = 1.
func TestVsatAboveCriticalTemperature(t *testing.T) {
	const (
		vc  = 255.0
		zc  = 0.274
		tol = 1e-12
	)

	for _, offset := range []float64{0.05, 0.2, 0.5} {
		below, err := liquids.Vsat(vc, zc, 1-offset)
		if err != nil {
			t.Fatalf("Vsat returned an unexpected error: %v", err)
		}

		above, err := liquids.Vsat(vc, zc, 1+offset)
		if err != nil {
			t.Fatalf("Vsat returned an unexpected error: %v", err)
		}

		if math.IsNaN(above) || above <= 0 {
			t.Errorf("at Tr = %.2f the volume is %v; want a finite positive value", 1+offset, above)
		}

		if math.Abs(above-below) > tol {
			t.Errorf(
				"the correlation is not symmetric about Tr = 1: %.9f at %.2f against %.9f at %.2f",
				above, 1+offset, below, 1-offset,
			)
		}
	}
}

// TestVsatInvalidInput checks the guards on the critical properties and
// the reduced temperature.
func TestVsatInvalidInput(t *testing.T) {
	testCases := []struct {
		name       string
		vc, zc, tr float64
	}{
		{"zero critical volume", 0, 0.274, 0.7},
		{"negative critical volume", -255, 0.274, 0.7},
		{"zero critical compressibility", 255, 0, 0.7},
		{"negative critical compressibility", 255, -0.274, 0.7},
		{"zero reduced temperature", 255, 0.274, 0},
		{"negative reduced temperature", 255, 0.274, -0.7},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := liquids.Vsat(tc.vc, tc.zc, tc.tr); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestReducedDensityIsPositive checks that the Lydersen chart returns a
// sensible density across the region it covers.
//
// Reduced density is the density divided by its critical value, so a
// compressed liquid sits above one and a gas below it. The values come
// from a digitised chart rather than an equation, so the check is on
// their being physically plausible rather than on any exact figure.
func TestReducedDensityIsPositive(t *testing.T) {
	var evaluated int

	for _, tr := range []float64{0.5, 0.6, 0.7, 0.8, 0.9, 1.0} {
		for _, pr := range []float64{0.5, 1, 2, 5, 10} {
			got, err := liquids.ReducedDensity(tr, pr)
			if err != nil {
				// Not every corner of the chart is populated.
				continue
			}

			evaluated++

			if got <= 0 || math.IsNaN(got) {
				t.Errorf("at Tr = %.1f and Pr = %.1f the reduced density is %v; want a positive value", tr, pr, got)
			}

			if got > 4 {
				t.Errorf("at Tr = %.1f and Pr = %.1f the reduced density is %.4f, which is implausibly high", tr, pr, got)
			}
		}
	}

	if evaluated == 0 {
		t.Fatal("no point of the chart could be evaluated")
	}
}

// TestReducedDensityFallsWithTemperature checks the direction of the
// chart: a liquid expands as it warms, so its reduced density falls
// along an isobar.
func TestReducedDensityFallsWithTemperature(t *testing.T) {
	const pr = 2.0

	var (
		previous  float64
		haveFirst bool
	)

	for _, tr := range []float64{0.5, 0.6, 0.7, 0.8, 0.9} {
		got, err := liquids.ReducedDensity(tr, pr)
		if err != nil {
			continue
		}

		if haveFirst && got >= previous {
			t.Errorf(
				"at Tr = %.1f the reduced density is %.4f; it should fall below the previous %.4f",
				tr, got, previous,
			)
		}

		previous, haveFirst = got, true
	}

	if !haveFirst {
		t.Fatal("no point of the isobar could be evaluated")
	}
}

// TestReducedDensityRisesWithPressure checks the other direction: a
// liquid is compressed by pressure, so its reduced density rises along
// an isotherm.
func TestReducedDensityRisesWithPressure(t *testing.T) {
	const tr = 0.7

	var (
		previous  float64
		haveFirst bool
	)

	for _, pr := range []float64{1, 2, 5, 10, 20} {
		got, err := liquids.ReducedDensity(tr, pr)
		if err != nil {
			continue
		}

		if haveFirst && got < previous {
			t.Errorf(
				"at Pr = %.1f the reduced density is %.4f; it should not fall below the previous %.4f",
				pr, got, previous,
			)
		}

		previous, haveFirst = got, true
	}

	if !haveFirst {
		t.Fatal("no point of the isotherm could be evaluated")
	}
}

// TestReducedDensityOutOfRange checks that conditions the chart does not
// cover are refused.
//
// The Lydersen chart is a digitisation of a printed figure and carries
// no information beyond its axes, so extrapolating would invent a
// number.
func TestReducedDensityOutOfRange(t *testing.T) {
	testCases := []struct {
		name   string
		tr, pr float64
	}{
		{"reduced temperature far below the chart", 0.01, 1},
		{"reduced temperature far above the chart", 100, 1},
		{"reduced pressure far below the chart", 0.7, 1e-6},
		{"reduced pressure far above the chart", 0.7, 1e6},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := liquids.ReducedDensity(tc.tr, tc.pr); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestReducedDensityInterpolatesBetweenIsotherms checks that a reduced
// temperature falling between two tabulated isotherms gives a value
// between theirs.
//
// The chart stores discrete isotherms and interpolates between them, so
// an intermediate temperature must not produce a result outside the pair
// that brackets it.
func TestReducedDensityInterpolatesBetweenIsotherms(t *testing.T) {
	const pr = 2.0

	lower, errLower := liquids.ReducedDensity(0.7, pr)
	upper, errUpper := liquids.ReducedDensity(0.8, pr)

	if errLower != nil || errUpper != nil {
		t.Skipf("the bracketing isotherms are unavailable at Pr = %g", pr)
	}

	middle, err := liquids.ReducedDensity(0.75, pr)
	if err != nil {
		t.Skipf("the intermediate isotherm is unavailable: %v", err)
	}

	low, high := math.Min(lower, upper), math.Max(lower, upper)

	if middle < low || middle > high {
		t.Errorf(
			"at Tr = 0.75 the reduced density is %.4f, outside the bracketing values %.4f and %.4f",
			middle, lower, upper,
		)
	}
}
