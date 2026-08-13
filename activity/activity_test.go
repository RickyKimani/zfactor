package activity_test

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor/activity"
	"github.com/rickykimani/zfactor/activity/margules"
	"github.com/rickykimani/zfactor/activity/nrtl"
	vanlaar "github.com/rickykimani/zfactor/activity/van-laar"
	"github.com/rickykimani/zfactor/activity/wilson"
)

// binaryModels returns one configured instance of every binary-capable
// model, so the shared invariants can be checked against all of them.
//
// The parameters are chosen to be strongly non-ideal and asymmetric:
// symmetric parameters make several formulations coincide and would
// hide index errors.
func binaryModels() map[string]activity.Model {
	return map[string]activity.Model{
		"Margules": margules.Margules{
			A12: 1.0, A21: 1.5,
			X: []float64{0.5, 0.5}, T: 318.15,
		},
		"VanLaar": vanlaar.VanLaar{
			A12: 1.0, A21: 1.5,
			X: []float64{0.5, 0.5}, T: 318.15,
		},
		"Wilson": wilson.Wilson{
			T: 318.15,
			X: []float64{0.5, 0.5},
			V: []float64{40.73, 79.84},
			Interaction: [][]float64{
				{0, 1000},
				{2000, 0},
			},
		},
		"NRTL": nrtl.NRTL{
			T: 318.15,
			X: []float64{0.5, 0.5},
			Alpha: [][]float64{
				{0, 0.3},
				{0.3, 0},
			},
			Tau: nrtl.ConstantTau{TauMatrix: [][]float64{
				{0, 0.5},
				{0.8, 0},
			}},
		},
	}
}

// TestPureComponentLimit checks that the activity coefficient of a
// component approaches unity as its mole fraction approaches one.
//
// A pure liquid is its own reference state, so this holds for every
// activity-coefficient model regardless of its parameters. It is the
// cheapest check that catches a misplaced composition index.
func TestPureComponentLimit(t *testing.T) {
	const (
		eps = 1e-10
		tol = 1e-9
	)

	for name, model := range binaryModels() {
		t.Run(name, func(t *testing.T) {
			for i := range 2 {
				x := []float64{eps, 1 - eps}
				if i == 0 {
					x = []float64{1 - eps, eps}
				}

				gamma, err := model.WithComposition(x).Activity()
				if err != nil {
					t.Fatalf("Activity returned an unexpected error: %v", err)
				}

				if math.Abs(gamma[i]-1) > tol {
					t.Errorf(
						"component %d at mole fraction %g: gamma = %.12f; want 1",
						i, 1-eps, gamma[i],
					)
				}
			}
		})
	}
}

// TestInfiniteDilutionMatchesLimit checks each model's closed-form
// infinite-dilution coefficients against its own Activity method
// evaluated at vanishing mole fraction.
//
// The two must agree: the closed form is the analytical limit of the
// same expression. They are nonetheless derived independently in the
// code, so a mistake in either is invisible until they are compared.
// activity.InfiniteDilution prefers the closed form whenever a model
// provides one, so an error there is returned to callers in place of
// the correct value rather than merely being unused.
func TestInfiniteDilutionMatchesLimit(t *testing.T) {
	const (
		eps = 1e-11
		tol = 1e-6
	)

	for name, model := range binaryModels() {
		t.Run(name, func(t *testing.T) {
			closed, ok := model.(activity.BinaryInfiniteDilutioner)
			if !ok {
				t.Skip("model provides no closed-form infinite dilution")
			}

			want, err := closed.BinaryInfiniteDilution()
			if err != nil {
				t.Fatalf("BinaryInfiniteDilution returned an unexpected error: %v", err)
			}

			for i := range 2 {
				x := []float64{1 - eps, eps}
				if i == 0 {
					x = []float64{eps, 1 - eps}
				}

				gamma, err := model.WithComposition(x).Activity()
				if err != nil {
					t.Fatalf("Activity returned an unexpected error: %v", err)
				}

				if rel := math.Abs(want[i]-gamma[i]) / gamma[i]; rel > tol {
					t.Errorf(
						"component %d: closed form gives gamma_inf = %.9f but the dilute limit is %.9f (%.2f%% apart)",
						i, want[i], gamma[i], 100*rel,
					)
				}
			}
		})
	}
}

// TestGibbsDuhem checks that each model satisfies the Gibbs-Duhem
// equation, which at constant temperature and pressure requires
//
//	x₁ d(ln γ₁)/dx₁ + x₂ d(ln γ₂)/dx₁ = 0.
//
// Any expression derived from a genuine excess Gibbs energy satisfies
// this identically, so a violation means the activity coefficients are
// not thermodynamically consistent — the signature of a wrong exponent,
// a dropped term or a transposed index.
//
// The derivatives are evaluated by central differences, so the residual
// is limited by the step size rather than by the models themselves.
func TestGibbsDuhem(t *testing.T) {
	const (
		h   = 1e-6
		tol = 1e-6
	)

	for name, model := range binaryModels() {
		t.Run(name, func(t *testing.T) {
			logGamma := func(x1 float64) (float64, float64) {
				t.Helper()

				gamma, err := model.
					WithComposition([]float64{x1, 1 - x1}).
					Activity()
				if err != nil {
					t.Fatalf("Activity returned an unexpected error: %v", err)
				}

				return math.Log(gamma[0]), math.Log(gamma[1])
			}

			for _, x1 := range []float64{0.2, 0.35, 0.5, 0.65, 0.8} {
				upper1, upper2 := logGamma(x1 + h)
				lower1, lower2 := logGamma(x1 - h)

				d1 := (upper1 - lower1) / (2 * h)
				d2 := (upper2 - lower2) / (2 * h)

				if residual := math.Abs(x1*d1 + (1-x1)*d2); residual > tol {
					t.Errorf(
						"at x1 = %.2f the Gibbs-Duhem residual is %.3e; want below %.0e",
						x1, residual, tol,
					)
				}
			}
		})
	}
}

// TestIdealParametersGiveUnitCoefficients checks that parameters
// describing an ideal solution produce unit activity coefficients at
// every composition, so that the models reduce to Raoult's law.
func TestIdealParametersGiveUnitCoefficients(t *testing.T) {
	const tol = 1e-12

	ideal := map[string]activity.Model{
		"Margules": margules.Margules{A12: 0, A21: 0, T: 318.15},
		"Wilson": wilson.Wilson{
			T: 318.15,
			V: []float64{50, 50},
			Interaction: [][]float64{
				{0, 0},
				{0, 0},
			},
		},
		"NRTL": nrtl.NRTL{
			T: 318.15,
			Alpha: [][]float64{
				{0, 0.3},
				{0.3, 0},
			},
			Tau: nrtl.ConstantTau{TauMatrix: [][]float64{
				{0, 0},
				{0, 0},
			}},
		},
	}

	// Van Laar is excluded: its parameters appear in denominators, so
	// zero values are rejected rather than describing an ideal solution.

	for name, model := range ideal {
		t.Run(name, func(t *testing.T) {
			for _, x1 := range []float64{0.1, 0.25, 0.5, 0.75, 0.9} {
				gamma, err := model.
					WithComposition([]float64{x1, 1 - x1}).
					Activity()
				if err != nil {
					t.Fatalf("Activity returned an unexpected error: %v", err)
				}

				for i, g := range gamma {
					if math.Abs(g-1) > tol {
						t.Errorf(
							"at x1 = %.2f: gamma%d = %.15f; want 1",
							x1, i+1, g,
						)
					}
				}
			}
		})
	}
}

// TestModelImmutability checks the contract documented on activity.Model:
// Composition must return a copy, and the With methods must leave the
// receiver untouched.
//
// The models are value types whose composition is a slice, so returning
// it directly would let a caller reach through and alter the model's
// state. The VLE solvers rely on this: they call WithComposition
// repeatedly while iterating and would corrupt their own inputs.
func TestModelImmutability(t *testing.T) {
	for name, model := range binaryModels() {
		t.Run(name, func(t *testing.T) {
			original := model.Composition()

			t.Run("Composition returns a copy", func(t *testing.T) {
				got := model.Composition()
				got[0] = 99

				if after := model.Composition(); after[0] != original[0] {
					t.Errorf(
						"mutating the returned slice changed the model: composition is now %v, was %v",
						after, original,
					)
				}
			})

			t.Run("WithComposition leaves the receiver unchanged", func(t *testing.T) {
				model.WithComposition([]float64{0.1, 0.9})

				if after := model.Composition(); after[0] != original[0] {
					t.Errorf(
						"receiver composition changed to %v; want %v",
						after, original,
					)
				}
			})

			t.Run("WithTemperature leaves the receiver unchanged", func(t *testing.T) {
				before := model.Temperature()
				model.WithTemperature(before + 100)

				if after := model.Temperature(); after != before {
					t.Errorf("receiver temperature changed to %g; want %g", after, before)
				}
			})

			t.Run("WithComposition copies its argument", func(t *testing.T) {
				x := []float64{0.3, 0.7}
				derived := model.WithComposition(x)
				x[0] = 99

				if got := derived.Composition(); got[0] != 0.3 {
					t.Errorf(
						"mutating the supplied slice changed the derived model: composition is now %v",
						got,
					)
				}
			})
		})
	}
}

// TestInfiniteDilutionGeneric checks the numerical routine in the root
// package against the closed forms the models provide.
//
// activity.InfiniteDilution dispatches to BinaryInfiniteDilution for
// binary mixtures that implement it, so it is also exercised here with
// a ternary mixture, where the numerical dilution path is the only one
// available.
func TestInfiniteDilutionGeneric(t *testing.T) {
	t.Run("binary dispatches to the closed form", func(t *testing.T) {
		const tol = 1e-9

		for name, model := range binaryModels() {
			closed, ok := model.(activity.BinaryInfiniteDilutioner)
			if !ok {
				continue
			}

			want, err := closed.BinaryInfiniteDilution()
			if err != nil {
				t.Fatalf("%s: BinaryInfiniteDilution returned an unexpected error: %v", name, err)
			}

			got, err := activity.InfiniteDilution(model)
			if err != nil {
				t.Fatalf("%s: InfiniteDilution returned an unexpected error: %v", name, err)
			}

			for i := range want {
				if math.Abs(got[i]-want[i]) > tol {
					t.Errorf("%s: component %d gamma_inf = %.9f; want %.9f", name, i, got[i], want[i])
				}
			}
		}
	})

	t.Run("ternary uses numerical dilution", func(t *testing.T) {
		model := wilson.Wilson{
			T: 318.15,
			X: []float64{1.0 / 3, 1.0 / 3, 1.0 / 3},
			V: []float64{40.73, 79.84, 60.0},
			Interaction: [][]float64{
				{0, 1000, 500},
				{2000, 0, 800},
				{300, 900, 0},
			},
		}

		got, err := activity.InfiniteDilution(model)
		if err != nil {
			t.Fatalf("InfiniteDilution returned an unexpected error: %v", err)
		}

		if len(got) != 3 {
			t.Fatalf("got %d coefficients; want 3", len(got))
		}

		for i, g := range got {
			if g <= 0 || math.IsNaN(g) || math.IsInf(g, 0) {
				t.Errorf("component %d: gamma_inf = %g; want a finite positive value", i, g)
			}
		}
	})

	t.Run("empty composition is rejected", func(t *testing.T) {
		model := wilson.Wilson{T: 318.15}

		if _, err := activity.InfiniteDilution(model); err == nil {
			t.Error("expected an error; got nil")
		}
	})
}

// TestCompositionValidation checks that malformed compositions are
// rejected by every model rather than producing a plausible-looking
// number.
func TestCompositionValidation(t *testing.T) {
	testCases := []struct {
		name string
		x    []float64
	}{
		{"does not sum to one", []float64{0.3, 0.3}},
		{"negative mole fraction", []float64{-0.2, 1.2}},
		{"wrong number of components", []float64{0.3, 0.3, 0.4}},
	}

	for name, model := range binaryModels() {
		t.Run(name, func(t *testing.T) {
			for _, tc := range testCases {
				t.Run(tc.name, func(t *testing.T) {
					if _, err := model.WithComposition(tc.x).Activity(); err == nil {
						t.Error("expected an error; got nil")
					}
				})
			}
		})
	}
}

// TestExtendedTau checks the temperature-dependent NRTL correlation
//
//	τij = aij + bij/T + cij ln(T) + dij T
//
// against a hand-evaluated value, and confirms that a non-positive
// temperature is rejected.
func TestExtendedTau(t *testing.T) {
	const T = 350.0

	correlation := nrtl.ExtendedTau{
		A: [][]float64{{0, 1.5}, {-0.5, 0}},
		B: [][]float64{{0, 200}, {300, 0}},
		C: [][]float64{{0, 0.1}, {0.2, 0}},
		D: [][]float64{{0, 0.001}, {0.002, 0}},
	}

	got, err := correlation.Tau(T)
	if err != nil {
		t.Fatalf("Tau returned an unexpected error: %v", err)
	}

	want := 1.5 + 200/T + 0.1*math.Log(T) + 0.001*T

	if math.Abs(got[0][1]-want) > 1e-12 {
		t.Errorf("tau12 = %.12f; want %.12f", got[0][1], want)
	}

	// The diagonal carries no parameters and must vanish.
	if got[0][0] != 0 {
		t.Errorf("tau11 = %g; want 0", got[0][0])
	}

	if _, err := correlation.Tau(0); err == nil {
		t.Error("expected an error for a non-positive temperature; got nil")
	}
}

// TestConstantTau checks that the constant correlation returns its
// matrix unchanged and ignores the temperature.
func TestConstantTau(t *testing.T) {
	want := [][]float64{{0, 0.5}, {0.8, 0}}

	correlation := nrtl.ConstantTau{TauMatrix: want}

	for _, temperature := range []float64{250, 318.15, 400} {
		got, err := correlation.Tau(temperature)
		if err != nil {
			t.Fatalf("Tau returned an unexpected error: %v", err)
		}

		for i := range want {
			for j := range want[i] {
				if got[i][j] != want[i][j] {
					t.Errorf(
						"at T = %g: tau[%d][%d] = %g; want %g",
						temperature, i, j, got[i][j], want[i][j],
					)
				}
			}
		}
	}
}

// TestMargulesExample13_1 checks the Margules coefficients against
// Example 13.1 of Smith, Van Ness & Abbott, for the
// methanol(1)/methyl acetate(2) system at 318.15 K.
//
// The example uses the symmetric form with
//
//	A = 2.771 - 0.00523 T,
//
// giving A = 1.10728 at this temperature. Part (a) evaluates the
// mixture at x1 = 0.25, where the published bubble pressure of
// 73.50 kPa follows from these coefficients and the saturation
// pressures of 44.51 and 65.64 kPa.
func TestMargulesExample13_1(t *testing.T) {
	const (
		a  = 2.771 - 0.00523*318.15
		x1 = 0.25

		wantGamma1 = 1.8641
		wantGamma2 = 1.0717

		tol = 1e-4
	)

	model := margules.Margules{A12: a, A21: a, T: 318.15}

	gamma, err := model.
		WithComposition([]float64{x1, 1 - x1}).
		Activity()
	if err != nil {
		t.Fatalf("Activity returned an unexpected error: %v", err)
	}

	if math.Abs(gamma[0]-wantGamma1) > tol {
		t.Errorf("gamma1 = %.6f; want %.4f", gamma[0], wantGamma1)
	}

	if math.Abs(gamma[1]-wantGamma2) > tol {
		t.Errorf("gamma2 = %.6f; want %.4f", gamma[1], wantGamma2)
	}

	// The published bubble pressure follows from the same coefficients.
	const (
		psat1    = 44.51
		psat2    = 65.64
		wantP    = 73.50
		pressTol = 0.01
	)

	p := x1*gamma[0]*psat1 + (1-x1)*gamma[1]*psat2

	if math.Abs(p-wantP) > pressTol {
		t.Errorf("bubble pressure = %.4f kPa; want %.2f kPa", p, wantP)
	}
}

// TestWilsonValidation checks the input validation specific to the
// Wilson model, whose parameters are matrices and molar volumes.
func TestWilsonValidation(t *testing.T) {
	base := wilson.Wilson{
		T: 318.15,
		X: []float64{0.5, 0.5},
		V: []float64{40.73, 79.84},
		Interaction: [][]float64{
			{0, 1000},
			{2000, 0},
		},
	}

	t.Run("non-positive temperature", func(t *testing.T) {
		model := base
		model.T = 0

		if _, err := model.Activity(); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("non-positive molar volume", func(t *testing.T) {
		model := base
		model.V = []float64{40.73, 0}

		if _, err := model.Activity(); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("wrong number of molar volumes", func(t *testing.T) {
		model := base
		model.V = []float64{40.73}

		if _, err := model.Activity(); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("non-square interaction matrix", func(t *testing.T) {
		model := base
		model.Interaction = [][]float64{{0, 1000}}

		if _, err := model.Activity(); err == nil {
			t.Error("expected an error; got nil")
		}
	})
}

// TestNRTLMissingTau checks that an NRTL model without a configured
// correlation is rejected rather than treated as ideal.
func TestNRTLMissingTau(t *testing.T) {
	model := nrtl.NRTL{
		T: 318.15,
		X: []float64{0.5, 0.5},
		Alpha: [][]float64{
			{0, 0.3},
			{0.3, 0},
		},
	}

	if _, err := model.Activity(); err == nil {
		t.Error("expected an error; got nil")
	}
}
