package modified_raoult

import (
	"fmt"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/activity/margules"
	"github.com/rickykimani/zfactor/activity/nrtl"
	vanlaar "github.com/rickykimani/zfactor/activity/van-laar"
	"github.com/rickykimani/zfactor/activity/wilson"
	"github.com/rickykimani/zfactor/vle/raoult"
)

// exampleTemperature is 318.15 K expressed in the °C this package works
// in, the condition of the isothermal parts of Example 13.1.
const exampleTemperature = 318.15 - 273.15

// exampleActivity returns the symmetric Margules model of Example 13.1
// evaluated at 318.15 K.
//
// The example's parameter varies with temperature, as A = 2.771 -
// 0.00523 T. Margules here holds its parameters fixed, so only the
// isothermal parts of the example can be reproduced; see
// TestIsobaricFlashesAreSelfConsistent for how the temperature solvers
// are checked instead.
func exampleActivity() ActivityModel {
	a := margulesA(318.15)
	return Margules{A12: a, A21: a}
}

// TestExample13_1IsothermalFlashes reproduces the two fixed-temperature
// parts of Example 13.1 of Smith, Van Ness & Abbott, for the
// methanol(1)/methyl acetate(2) system at 318.15 K.
//
// Methyl acetate is the more volatile component, so the vapor is richer
// in it and correspondingly the liquid in methanol. Part (b) is checked
// against the equilibrium relation itself rather than a printed value,
// since only the pressure is quoted in the text.
func TestExample13_1IsothermalFlashes(t *testing.T) {
	t.Run("bubble pressure", func(t *testing.T) {
		const (
			wantP  = 73.50
			wantY1 = 0.282
		)

		got, err := BubbleP(MixtureInput{
			T:            exampleTemperature,
			Compositions: []float64{0.25, 0.75},
			Antoine:      exampleModels(),
			Activity:     exampleActivity(),
		})
		if err != nil {
			t.Fatalf("BubbleP returned an unexpected error: %v", err)
		}

		if math.Abs(got.P-wantP) > 0.01 {
			t.Errorf("bubble pressure = %.4f kPa; want %.2f", got.P, wantP)
		}

		if math.Abs(got.Y[0]-wantY1) > 5e-4 {
			t.Errorf("vapor composition y1 = %.4f; want %.3f", got.Y[0], wantY1)
		}
	})

	t.Run("dew pressure", func(t *testing.T) {
		const wantP = 62.89

		got, err := DewP(MixtureInput{
			T:            exampleTemperature,
			Compositions: []float64{0.60, 0.40},
			Antoine:      exampleModels(),
			Activity:     exampleActivity(),
		})
		if err != nil {
			t.Fatalf("DewP returned an unexpected error: %v", err)
		}

		if math.Abs(got.P-wantP) > 0.01 {
			t.Errorf("dew pressure = %.4f kPa; want %.2f", got.P, wantP)
		}

		// Methanol is the less volatile component, so the liquid holds
		// more of it than the vapor does.
		if got.X[0] <= 0.60 {
			t.Errorf("liquid composition x1 = %.4f; want more methanol than the vapor's 0.60", got.X[0])
		}

		assertEquilibrium(t, got.P, got.X, []float64{0.60, 0.40}, exampleTemperature, exampleActivity())
	})
}

// assertEquilibrium checks the relation the package solves, that the
// partial pressure of each component matches on both sides of the phase
// boundary:
//
//	yi P = xi gamma_i Pi_sat.
//
// It is the definition of the answer rather than a value read from
// anywhere, so it holds at every converged state.
//
// The comparison is relative. Where the temperature was solved for, the
// iteration converges on the temperature itself, and the residual left
// in the partial pressures is that tolerance multiplied by dPsat/dT,
// which for these substances is a few kPa per kelvin.
func assertEquilibrium(
	t *testing.T,
	pressure float64,
	x, y []float64,
	temperature float64,
	model ActivityModel,
) {
	t.Helper()

	const relTol = 1e-5

	gamma, err := activityCoefficients(model, temperature, x)
	if err != nil {
		t.Fatalf("activityCoefficients returned an unexpected error: %v", err)
	}

	psat, err := MixtureInput{T: temperature, Antoine: exampleModels()}.PSat()
	if err != nil {
		t.Fatalf("PSat returned an unexpected error: %v", err)
	}

	for i := range x {
		vapor := y[i] * pressure
		liquid := x[i] * gamma[i] * psat[i]

		scale := math.Max(math.Abs(vapor), 1)

		if rel := math.Abs(vapor-liquid) / scale; rel > relTol {
			t.Errorf(
				"component %d: y*P = %.9f but x*gamma*Psat = %.9f (%.2e apart, relative)",
				i, vapor, liquid, rel,
			)
		}
	}
}

// TestIdealActivityReducesToRaoult checks that the modified law collapses
// onto the simple one when the solution is ideal.
//
// Margules parameters of zero give unit activity coefficients at every
// composition, at which point yi P = xi Pi_sat, which is Raoult's law.
// The two packages implement that separately, so agreement to rounding
// is evidence that the activity coefficients enter in the right places.
func TestIdealActivityReducesToRaoult(t *testing.T) {
	const tol = 1e-9

	ideal := Margules{A12: 0, A21: 0}
	models := exampleModels()

	for _, z1 := range []float64{0.1, 0.25, 0.5, 0.75, 0.9} {
		composition := []float64{z1, 1 - z1}

		t.Run("bubble pressure", func(t *testing.T) {
			modified, err := BubbleP(MixtureInput{
				T: exampleTemperature, Compositions: composition,
				Antoine: models, Activity: ideal,
			})
			if err != nil {
				t.Fatalf("BubbleP returned an unexpected error: %v", err)
			}

			simple, err := raoult.BubbleP(raoult.MixtureInput{
				T: exampleTemperature, Compositions: composition, Antoine: models,
			})
			if err != nil {
				t.Fatalf("raoult.BubbleP returned an unexpected error: %v", err)
			}

			if math.Abs(modified.P-simple.P) > tol {
				t.Errorf("at z1 = %.2f: %.10f kPa against Raoult's %.10f", z1, modified.P, simple.P)
			}

			if math.Abs(modified.Y[0]-simple.Y[0]) > tol {
				t.Errorf("at z1 = %.2f: y1 = %.10f against Raoult's %.10f", z1, modified.Y[0], simple.Y[0])
			}
		})

		t.Run("dew pressure", func(t *testing.T) {
			modified, err := DewP(MixtureInput{
				T: exampleTemperature, Compositions: composition,
				Antoine: models, Activity: ideal,
			})
			if err != nil {
				t.Fatalf("DewP returned an unexpected error: %v", err)
			}

			simple, err := raoult.DewP(raoult.MixtureInput{
				T: exampleTemperature, Compositions: composition, Antoine: models,
			})
			if err != nil {
				t.Fatalf("raoult.DewP returned an unexpected error: %v", err)
			}

			if math.Abs(modified.P-simple.P) > 1e-6 {
				t.Errorf("at z1 = %.2f: %.10f kPa against Raoult's %.10f", z1, modified.P, simple.P)
			}

			if math.Abs(modified.X[0]-simple.X[0]) > 1e-6 {
				t.Errorf("at z1 = %.2f: x1 = %.10f against Raoult's %.10f", z1, modified.X[0], simple.X[0])
			}
		})
	}
}

// TestBubbleDewPressureRoundTrip checks that the two isothermal flashes
// describe the same phase boundary from opposite sides.
//
// Feeding the vapor composition from a bubble calculation into a dew
// calculation at the same temperature must return the original pressure
// and liquid composition, since both name the same equilibrium state.
// Unlike Raoult's law this is not an algebraic identity — the dew
// calculation iterates on the activity coefficients — so the agreement
// is bounded by the solver tolerance.
func TestBubbleDewPressureRoundTrip(t *testing.T) {
	const tol = 1e-6

	for _, x1 := range []float64{0.1, 0.25, 0.5, 0.75, 0.9} {
		bubble, err := BubbleP(MixtureInput{
			T: exampleTemperature, Compositions: []float64{x1, 1 - x1},
			Antoine: exampleModels(), Activity: exampleActivity(),
		})
		if err != nil {
			t.Fatalf("BubbleP returned an unexpected error: %v", err)
		}

		dew, err := DewP(MixtureInput{
			T: exampleTemperature, Compositions: bubble.Y,
			Antoine: exampleModels(), Activity: exampleActivity(),
		})
		if err != nil {
			t.Fatalf("DewP returned an unexpected error: %v", err)
		}

		if math.Abs(dew.P-bubble.P) > tol {
			t.Errorf(
				"at x1 = %.2f: the dew pressure %.9f does not recover the bubble pressure %.9f",
				x1, dew.P, bubble.P,
			)
		}

		if math.Abs(dew.X[0]-x1) > tol {
			t.Errorf("at x1 = %.2f: the recovered liquid composition is %.9f", x1, dew.X[0])
		}
	}
}

// TestIsobaricFlashesAreSelfConsistent checks the temperature solvers
// against the pressure solvers.
//
// The isobaric parts of Example 13.1 cannot be reproduced directly,
// because the example's Margules parameter varies with temperature while
// this package's model holds it fixed. Solving at a fixed parameter is
// still a well-posed problem, so the temperature returned must satisfy
// the corresponding pressure calculation: a bubble temperature fed back
// into BubbleP must return the pressure it was found for.
func TestIsobaricFlashesAreSelfConsistent(t *testing.T) {
	const (
		pressure = 101.33
		tol      = 1e-6
	)

	for _, z1 := range []float64{0.15, 0.4, 0.6, 0.85} {
		composition := []float64{z1, 1 - z1}

		t.Run("bubble", func(t *testing.T) {
			bubbleT, err := BubbleT(MixtureInput{
				P: pressure, Compositions: composition,
				Antoine: exampleModels(), Activity: exampleActivity(),
			})
			if err != nil {
				t.Fatalf("BubbleT returned an unexpected error: %v", err)
			}

			back, err := BubbleP(MixtureInput{
				T: bubbleT.T, Compositions: composition,
				Antoine: exampleModels(), Activity: exampleActivity(),
			})
			if err != nil {
				t.Fatalf("BubbleP returned an unexpected error: %v", err)
			}

			if math.Abs(back.P-pressure) > tol {
				t.Errorf(
					"at x1 = %.2f: the bubble temperature %.6f °C reproduces %.9f kPa",
					z1, bubbleT.T, back.P,
				)
			}

			assertEquilibrium(t, pressure, composition, bubbleT.Y, bubbleT.T, exampleActivity())
		})

		t.Run("dew", func(t *testing.T) {
			dewT, err := DewT(MixtureInput{
				P: pressure, Compositions: composition,
				Antoine: exampleModels(), Activity: exampleActivity(),
			})
			if err != nil {
				t.Fatalf("DewT returned an unexpected error: %v", err)
			}

			back, err := DewP(MixtureInput{
				T: dewT.T, Compositions: composition,
				Antoine: exampleModels(), Activity: exampleActivity(),
			})
			if err != nil {
				t.Fatalf("DewP returned an unexpected error: %v", err)
			}

			if math.Abs(back.P-pressure) > tol {
				t.Errorf(
					"at y1 = %.2f: the dew temperature %.6f °C reproduces %.9f kPa",
					z1, dewT.T, back.P,
				)
			}

			assertEquilibrium(t, pressure, dewT.X, composition, dewT.T, exampleActivity())
		})
	}
}

// TestCompositionsAreNormalised checks that every composition the
// package reports sums to one.
func TestCompositionsAreNormalised(t *testing.T) {
	const tol = 1e-9

	composition := []float64{0.4, 0.6}

	total := func(v []float64) float64 {
		var sum float64
		for _, e := range v {
			sum += e
		}
		return sum
	}

	bubbleP, err := BubbleP(MixtureInput{
		T: exampleTemperature, Compositions: composition,
		Antoine: exampleModels(), Activity: exampleActivity(),
	})
	if err != nil {
		t.Fatalf("BubbleP returned an unexpected error: %v", err)
	}
	if math.Abs(total(bubbleP.Y)-1) > tol {
		t.Errorf("BubbleP vapor composition sums to %.12f", total(bubbleP.Y))
	}

	dewP, err := DewP(MixtureInput{
		T: exampleTemperature, Compositions: composition,
		Antoine: exampleModels(), Activity: exampleActivity(),
	})
	if err != nil {
		t.Fatalf("DewP returned an unexpected error: %v", err)
	}
	if math.Abs(total(dewP.X)-1) > tol {
		t.Errorf("DewP liquid composition sums to %.12f", total(dewP.X))
	}

	bubbleT, err := BubbleT(MixtureInput{
		P: 101.33, Compositions: composition,
		Antoine: exampleModels(), Activity: exampleActivity(),
	})
	if err != nil {
		t.Fatalf("BubbleT returned an unexpected error: %v", err)
	}
	if math.Abs(total(bubbleT.Y)-1) > tol {
		t.Errorf("BubbleT vapor composition sums to %.12f", total(bubbleT.Y))
	}

	dewT, err := DewT(MixtureInput{
		P: 101.33, Compositions: composition,
		Antoine: exampleModels(), Activity: exampleActivity(),
	})
	if err != nil {
		t.Fatalf("DewT returned an unexpected error: %v", err)
	}
	if math.Abs(total(dewT.X)-1) > tol {
		t.Errorf("DewT liquid composition sums to %.12f", total(dewT.X))
	}
}

// TestActivityModelWrappers checks that each wrapper constructs the
// underlying model from the activity package.
//
// The wrappers exist so that the VLE calculations depend on this
// package's own types rather than on the activity implementations, and
// each simply forwards its parameters.
func TestActivityModelWrappers(t *testing.T) {
	testCases := []struct {
		name    string
		wrapper ActivityModel
		want    any
	}{
		{"Margules", Margules{A12: 1, A21: 1.5}, margules.Margules{}},
		{"VanLaar", VanLaar{A12: 1, A21: 1.5}, vanlaar.VanLaar{}},
		{
			name: "Wilson",
			wrapper: Wilson{
				V:           []float64{40.73, 79.84},
				Interaction: [][]float64{{0, 1000}, {2000, 0}},
			},
			want: wilson.Wilson{},
		},
		{
			name: "NRTL",
			wrapper: NRTL{
				Alpha: [][]float64{{0, 0.3}, {0.3, 0}},
				Tau:   nrtl.ConstantTau{TauMatrix: [][]float64{{0, 0.5}, {0.8, 0}}},
			},
			want: nrtl.NRTL{},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got := tc.wrapper.Model()

			if gotType, wantType := fmt.Sprintf("%T", got), fmt.Sprintf("%T", tc.want); gotType != wantType {
				t.Fatalf("Model() returned %s; want %s", gotType, wantType)
			}

			// The constructed model must be usable at a real state.
			gamma, err := got.
				WithTemperature(318.15).
				WithComposition([]float64{0.4, 0.6}).
				Activity()
			if err != nil {
				t.Fatalf("Activity returned an unexpected error: %v", err)
			}

			for i, g := range gamma {
				if g <= 0 || math.IsNaN(g) {
					t.Errorf("gamma%d = %g; want a positive value", i+1, g)
				}
			}
		})
	}
}

// TestInvalidInput checks that malformed mixtures are rejected by every
// entry point.
func TestInvalidInput(t *testing.T) {
	testCases := []struct {
		name  string
		input MixtureInput
	}{
		{
			name: "composition does not sum to one",
			input: MixtureInput{
				T: exampleTemperature, P: 101.33,
				Compositions: []float64{0.3, 0.3},
				Antoine:      exampleModels(), Activity: exampleActivity(),
			},
		},
		{
			name: "negative mole fraction",
			input: MixtureInput{
				T: exampleTemperature, P: 101.33,
				Compositions: []float64{-0.2, 1.2},
				Antoine:      exampleModels(), Activity: exampleActivity(),
			},
		},
		{
			name: "no activity model",
			input: MixtureInput{
				T: exampleTemperature, P: 101.33,
				Compositions: []float64{0.4, 0.6},
				Antoine:      exampleModels(),
			},
		},
		{
			name: "more components than correlations",
			input: MixtureInput{
				T: exampleTemperature, P: 101.33,
				Compositions: []float64{0.3, 0.3, 0.4},
				Antoine:      exampleModels(), Activity: exampleActivity(),
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := BubbleP(tc.input); err == nil {
				t.Error("BubbleP: expected an error; got nil")
			}

			if _, err := DewP(tc.input); err == nil {
				t.Error("DewP: expected an error; got nil")
			}

			if _, err := BubbleT(tc.input); err == nil {
				t.Error("BubbleT: expected an error; got nil")
			}

			if _, err := DewT(tc.input); err == nil {
				t.Error("DewT: expected an error; got nil")
			}
		})
	}
}
