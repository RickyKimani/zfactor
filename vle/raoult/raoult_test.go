package raoult_test

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/vle/raoult"
)

// benzeneToluene returns the Antoine correlations for the benzene(1) /
// toluene(2) system, the standard worked example for Raoult's law.
//
// The tabulated coefficients are those of Table B.2, so results here are
// directly comparable with hand calculations made from the same table.
func benzeneToluene() []antoine.Model {
	return []antoine.Model{antoine.Benzene, antoine.Toluene}
}

// TestBenzeneTolueneReference reproduces three hand calculations for the
// benzene(1)/toluene(2) system, one for each kind of flash the package
// supports.
//
// Benzene is the more volatile component throughout, so the vapor is
// always richer in it than the liquid.
func TestBenzeneTolueneReference(t *testing.T) {
	const (
		compositionTol = 5e-4
		pressureTol    = 5e-3
		temperatureTol = 5e-3
	)

	t.Run("bubble pressure", func(t *testing.T) {
		// x1 = 0.33 at 100 °C.
		got, err := raoult.BubbleP(raoult.MixtureInput{
			T:            100,
			Compositions: []float64{0.33, 0.67},
			Antoine:      benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("BubbleP returned an unexpected error: %v", err)
		}

		if math.Abs(got.P-109.303) > pressureTol {
			t.Errorf("bubble pressure = %.4f kPa; want 109.303", got.P)
		}

		if math.Abs(got.Y[0]-0.545) > compositionTol {
			t.Errorf("vapor composition y1 = %.4f; want 0.545", got.Y[0])
		}
	})

	t.Run("dew pressure", func(t *testing.T) {
		// y1 = 0.33 at 100 °C.
		got, err := raoult.DewP(raoult.MixtureInput{
			T:            100,
			Compositions: []float64{0.33, 0.67},
			Antoine:      benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("DewP returned an unexpected error: %v", err)
		}

		if math.Abs(got.P-92.156) > pressureTol {
			t.Errorf("dew pressure = %.4f kPa; want 92.156", got.P)
		}

		if math.Abs(got.X[0]-0.169) > compositionTol {
			t.Errorf("liquid composition x1 = %.4f; want 0.169", got.X[0])
		}
	})

	t.Run("bubble temperature", func(t *testing.T) {
		// x1 = 0.33 at 120 kPa.
		got, err := raoult.BubbleT(raoult.MixtureInput{
			P:            120,
			Compositions: []float64{0.33, 0.67},
			Antoine:      benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("BubbleT returned an unexpected error: %v", err)
		}

		if math.Abs(got.T-103.307) > temperatureTol {
			t.Errorf("bubble temperature = %.4f °C; want 103.307", got.T)
		}

		if math.Abs(got.Y[0]-0.542) > compositionTol {
			t.Errorf("vapor composition y1 = %.4f; want 0.542", got.Y[0])
		}
	})
}

// TestBubbleDewPressureRoundTrip checks an exact identity of Raoult's law.
//
// Bubble pressure takes the composition-weighted arithmetic mean of the
// saturation pressures and dew pressure the harmonic mean. Feeding the
// vapor composition from a bubble calculation into a dew calculation at
// the same temperature must return both the original pressure and the
// original liquid composition, since
//
//	Σ yi/Pi_sat = Σ xi Pi_sat / (P Pi_sat) = 1/P.
//
// The relation is algebraic rather than iterative, so it holds to
// rounding for any composition and temperature.
func TestBubbleDewPressureRoundTrip(t *testing.T) {
	const tol = 1e-9

	for _, x1 := range []float64{0.05, 0.2, 0.33, 0.5, 0.75, 0.95} {
		bubble, err := raoult.BubbleP(raoult.MixtureInput{
			T:            100,
			Compositions: []float64{x1, 1 - x1},
			Antoine:      benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("BubbleP returned an unexpected error: %v", err)
		}

		dew, err := raoult.DewP(raoult.MixtureInput{
			T:            100,
			Compositions: bubble.Y,
			Antoine:      benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("DewP returned an unexpected error: %v", err)
		}

		if math.Abs(dew.P-bubble.P) > tol {
			t.Errorf(
				"at x1 = %.2f: dew pressure %.12f does not recover the bubble pressure %.12f",
				x1, dew.P, bubble.P,
			)
		}

		if math.Abs(dew.X[0]-x1) > tol {
			t.Errorf(
				"at x1 = %.2f: recovered liquid composition %.12f",
				x1, dew.X[0],
			)
		}
	}
}

// TestBubbleTemperatureRoundTrip checks the temperature and pressure
// solvers against each other.
//
// The bubble temperature at a given pressure must, when fed back into
// the bubble-pressure calculation, return that pressure. The two solve
// the same equation for different unknowns, so agreement is limited only
// by the tolerance of the iteration.
func TestBubbleTemperatureRoundTrip(t *testing.T) {
	const tol = 1e-6

	for _, x1 := range []float64{0.1, 0.33, 0.5, 0.9} {
		for _, pressure := range []float64{80, 101.325, 120} {
			bubbleT, err := raoult.BubbleT(raoult.MixtureInput{
				P:            pressure,
				Compositions: []float64{x1, 1 - x1},
				Antoine:      benzeneToluene(),
			})
			if err != nil {
				t.Fatalf("BubbleT returned an unexpected error: %v", err)
			}

			bubbleP, err := raoult.BubbleP(raoult.MixtureInput{
				T:            bubbleT.T,
				Compositions: []float64{x1, 1 - x1},
				Antoine:      benzeneToluene(),
			})
			if err != nil {
				t.Fatalf("BubbleP returned an unexpected error: %v", err)
			}

			if math.Abs(bubbleP.P-pressure) > tol {
				t.Errorf(
					"at x1 = %.2f and P = %g kPa: the bubble temperature %.6f °C reproduces %.9f kPa",
					x1, pressure, bubbleT.T, bubbleP.P,
				)
			}
		}
	}
}

// TestDewTemperatureRoundTrip applies the same check to the dew solvers.
func TestDewTemperatureRoundTrip(t *testing.T) {
	const tol = 1e-6

	for _, y1 := range []float64{0.1, 0.33, 0.5, 0.9} {
		dewT, err := raoult.DewT(raoult.MixtureInput{
			P:            101.325,
			Compositions: []float64{y1, 1 - y1},
			Antoine:      benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("DewT returned an unexpected error: %v", err)
		}

		dewP, err := raoult.DewP(raoult.MixtureInput{
			T:            dewT.T,
			Compositions: []float64{y1, 1 - y1},
			Antoine:      benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("DewP returned an unexpected error: %v", err)
		}

		if math.Abs(dewP.P-101.325) > tol {
			t.Errorf(
				"at y1 = %.2f: the dew temperature %.6f °C reproduces %.9f kPa",
				y1, dewT.T, dewP.P,
			)
		}
	}
}

// TestPureComponentLimit checks that a mixture of one component reduces
// to that component's saturation pressure.
//
// With a single species present, Raoult's law states P = Pi_sat and the
// vapor has the same composition as the liquid, whichever flash is used.
func TestPureComponentLimit(t *testing.T) {
	const tol = 1e-9

	testCases := []struct {
		name  string
		x     []float64
		model *antoine.Antoine
	}{
		{"pure benzene", []float64{1, 0}, antoine.Benzene},
		{"pure toluene", []float64{0, 1}, antoine.Toluene},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			want, err := tc.model.Pressure(100)
			if err != nil {
				t.Fatalf("Pressure returned an unexpected error: %v", err)
			}

			bubble, err := raoult.BubbleP(raoult.MixtureInput{
				T:            100,
				Compositions: tc.x,
				Antoine:      benzeneToluene(),
			})
			if err != nil {
				t.Fatalf("BubbleP returned an unexpected error: %v", err)
			}

			if math.Abs(bubble.P-want) > tol {
				t.Errorf("bubble pressure = %.9f kPa; want the pure saturation pressure %.9f", bubble.P, want)
			}

			for i := range tc.x {
				if math.Abs(bubble.Y[i]-tc.x[i]) > tol {
					t.Errorf("vapor composition = %v; want the liquid composition %v", bubble.Y, tc.x)
					break
				}
			}
		})
	}
}

// TestBubblePressureExceedsDewPressure checks the ordering of the two
// pressures at a fixed composition.
//
// The bubble pressure is a weighted arithmetic mean of the saturation
// pressures and the dew pressure the corresponding harmonic mean, so the
// first can never fall below the second. Equality would require every
// component to share a saturation pressure.
func TestBubblePressureExceedsDewPressure(t *testing.T) {
	for _, z1 := range []float64{0.1, 0.25, 0.5, 0.75, 0.9} {
		composition := []float64{z1, 1 - z1}

		bubble, err := raoult.BubbleP(raoult.MixtureInput{
			T: 100, Compositions: composition, Antoine: benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("BubbleP returned an unexpected error: %v", err)
		}

		dew, err := raoult.DewP(raoult.MixtureInput{
			T: 100, Compositions: composition, Antoine: benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("DewP returned an unexpected error: %v", err)
		}

		if bubble.P <= dew.P {
			t.Errorf(
				"at z1 = %.2f: bubble pressure %.6f must exceed dew pressure %.6f",
				z1, bubble.P, dew.P,
			)
		}
	}
}

// TestBubbleTemperatureBelowDewTemperature checks the matching ordering
// in temperature: at a given pressure a mixture begins to boil below the
// temperature at which it finishes condensing.
func TestBubbleTemperatureBelowDewTemperature(t *testing.T) {
	for _, z1 := range []float64{0.1, 0.33, 0.5, 0.9} {
		composition := []float64{z1, 1 - z1}

		bubble, err := raoult.BubbleT(raoult.MixtureInput{
			P: 101.325, Compositions: composition, Antoine: benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("BubbleT returned an unexpected error: %v", err)
		}

		dew, err := raoult.DewT(raoult.MixtureInput{
			P: 101.325, Compositions: composition, Antoine: benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("DewT returned an unexpected error: %v", err)
		}

		if bubble.T >= dew.T {
			t.Errorf(
				"at z1 = %.2f: bubble temperature %.6f must lie below dew temperature %.6f",
				z1, bubble.T, dew.T,
			)
		}
	}
}

// TestVaporIsRicherInTheVolatileComponent checks the direction of the
// separation that makes distillation work.
//
// Benzene has the higher saturation pressure at these conditions, so the
// vapor leaving a boiling liquid must carry more of it than the liquid
// holds, and the liquid condensing from a vapor must carry less.
func TestVaporIsRicherInTheVolatileComponent(t *testing.T) {
	for _, z1 := range []float64{0.1, 0.25, 0.5, 0.75, 0.9} {
		composition := []float64{z1, 1 - z1}

		bubble, err := raoult.BubbleP(raoult.MixtureInput{
			T: 100, Compositions: composition, Antoine: benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("BubbleP returned an unexpected error: %v", err)
		}

		if bubble.Y[0] <= z1 {
			t.Errorf(
				"at x1 = %.2f: vapor benzene fraction %.4f must exceed the liquid fraction",
				z1, bubble.Y[0],
			)
		}

		dew, err := raoult.DewP(raoult.MixtureInput{
			T: 100, Compositions: composition, Antoine: benzeneToluene(),
		})
		if err != nil {
			t.Fatalf("DewP returned an unexpected error: %v", err)
		}

		if dew.X[0] >= z1 {
			t.Errorf(
				"at y1 = %.2f: liquid benzene fraction %.4f must fall below the vapor fraction",
				z1, dew.X[0],
			)
		}
	}
}

// TestEquilibriumCompositionsSumToOne checks that every composition the
// package reports is normalised.
func TestEquilibriumCompositionsSumToOne(t *testing.T) {
	const tol = 1e-9

	composition := []float64{0.33, 0.67}

	sum := func(v []float64) float64 {
		var total float64
		for _, e := range v {
			total += e
		}
		return total
	}

	bubbleP, err := raoult.BubbleP(raoult.MixtureInput{
		T: 100, Compositions: composition, Antoine: benzeneToluene(),
	})
	if err != nil {
		t.Fatalf("BubbleP returned an unexpected error: %v", err)
	}
	if math.Abs(sum(bubbleP.Y)-1) > tol {
		t.Errorf("BubbleP vapor composition sums to %.12f; want 1", sum(bubbleP.Y))
	}

	dewP, err := raoult.DewP(raoult.MixtureInput{
		T: 100, Compositions: composition, Antoine: benzeneToluene(),
	})
	if err != nil {
		t.Fatalf("DewP returned an unexpected error: %v", err)
	}
	if math.Abs(sum(dewP.X)-1) > tol {
		t.Errorf("DewP liquid composition sums to %.12f; want 1", sum(dewP.X))
	}

	bubbleT, err := raoult.BubbleT(raoult.MixtureInput{
		P: 101.325, Compositions: composition, Antoine: benzeneToluene(),
	})
	if err != nil {
		t.Fatalf("BubbleT returned an unexpected error: %v", err)
	}
	if math.Abs(sum(bubbleT.Y)-1) > tol {
		t.Errorf("BubbleT vapor composition sums to %.12f; want 1", sum(bubbleT.Y))
	}

	dewT, err := raoult.DewT(raoult.MixtureInput{
		P: 101.325, Compositions: composition, Antoine: benzeneToluene(),
	})
	if err != nil {
		t.Fatalf("DewT returned an unexpected error: %v", err)
	}
	if math.Abs(sum(dewT.X)-1) > tol {
		t.Errorf("DewT liquid composition sums to %.12f; want 1", sum(dewT.X))
	}
}

// TestSuppliedSaturationPressures checks the input type that bypasses the
// Antoine correlations, so that a caller with measured saturation
// pressures can use the same solvers.
func TestSuppliedSaturationPressures(t *testing.T) {
	const tol = 1e-9

	// The benzene and toluene saturation pressures at 100 °C.
	psat := []float64{180.4528, 74.2597}
	x := []float64{0.33, 0.67}

	got, err := raoult.BubbleP(raoult.SaturationPressureInput{
		Compositions: x,
		PSats:        psat,
	})
	if err != nil {
		t.Fatalf("BubbleP returned an unexpected error: %v", err)
	}

	want := x[0]*psat[0] + x[1]*psat[1]

	if math.Abs(got.P-want) > tol {
		t.Errorf("bubble pressure = %.9f; want the weighted mean %.9f", got.P, want)
	}
}

// TestInvalidInput checks that malformed mixtures are rejected rather
// than producing a plausible-looking number.
func TestInvalidInput(t *testing.T) {
	testCases := []struct {
		name  string
		input raoult.MixtureInput
	}{
		{
			name: "composition does not sum to one",
			input: raoult.MixtureInput{
				T: 100, Compositions: []float64{0.3, 0.3}, Antoine: benzeneToluene(),
			},
		},
		{
			name: "negative mole fraction",
			input: raoult.MixtureInput{
				T: 100, Compositions: []float64{-0.2, 1.2}, Antoine: benzeneToluene(),
			},
		},
		{
			name: "more components than correlations",
			input: raoult.MixtureInput{
				T: 100, Compositions: []float64{0.3, 0.3, 0.4}, Antoine: benzeneToluene(),
			},
		},
		{
			name: "no components",
			input: raoult.MixtureInput{
				T: 100, Compositions: nil, Antoine: benzeneToluene(),
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := raoult.BubbleP(tc.input); err == nil {
				t.Error("BubbleP: expected an error; got nil")
			}

			if _, err := raoult.DewP(tc.input); err == nil {
				t.Error("DewP: expected an error; got nil")
			}
		})
	}
}
