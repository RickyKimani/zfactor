package raoult_test

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/vle"
	"github.com/rickykimani/zfactor/vle/raoult"
)

// tight converges far below the accuracy any assertion here needs, so a
// failure indicates a defect rather than a tolerance.
func tight() vle.SolverOptions {
	return vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 200}
}

// TestFlashExample13_8 reproduces Example 13.8 of Smith, Van Ness &
// Abbott: acetone(1)/acetonitrile(2)/nitromethane(3) at 80 °C and
// 110 kPa, with an overall composition of 0.45, 0.35 and 0.20.
//
// The saturation pressures are supplied directly, as the example does,
// rather than taken from a correlation.
func TestFlashExample13_8(t *testing.T) {
	got, err := raoult.FlashPT(raoult.SaturationPressureInput{
		P:            110,
		Compositions: []float64{0.45, 0.35, 0.20},
		PSats:        []float64{195.75, 97.84, 50.32},
		Options:      tight(),
	})
	if err != nil {
		t.Fatalf("FlashPT returned an unexpected error: %v", err)
	}

	if math.Abs(got.V-0.7364) > 5e-4 {
		t.Errorf("vapor = %.6f mol; want 0.7364", got.V)
	}

	if math.Abs(got.L-0.2636) > 5e-4 {
		t.Errorf("liquid = %.6f mol; want 0.2636", got.L)
	}

	if math.Abs(got.V+got.L-1) > 1e-12 {
		t.Errorf("the phases hold %.12f mol between them; want 1", got.V+got.L)
	}

	for i, want := range []float64{0.5087, 0.3389, 0.1524} {
		if math.Abs(got.Y[i]-want) > 5e-4 {
			t.Errorf("y%d = %.4f; want %.4f", i+1, got.Y[i], want)
		}
	}

	for i, want := range []float64{0.2859, 0.3810, 0.3331} {
		if math.Abs(got.X[i]-want) > 5e-4 {
			t.Errorf("x%d = %.4f; want %.4f", i+1, got.X[i], want)
		}
	}
}

// TestFlashMatchesTheBoundariesAtTheirPressures ties the flash to the
// bubble- and dew-point calculations.
//
// A feed at its bubble pressure has only begun to evaporate, so the vapor
// fraction is zero; at its dew pressure only the last drop of liquid
// remains, so the fraction is one. The three calculations are separate
// code paths solving the same equilibrium, and these two conditions are
// where they must agree.
//
// Under Raoult's law the boundaries follow from the feed composition
// alone, which is what makes the comparison exact rather than
// approximate.
func TestFlashMatchesTheBoundariesAtTheirPressures(t *testing.T) {
	const tol = 1e-6

	z := []float64{0.4, 0.35, 0.25}
	psat := []float64{195.75, 97.84, 50.32}

	bubble, err := raoult.BubbleP(raoult.SaturationPressureInput{
		Compositions: z, PSats: psat,
	})
	if err != nil {
		t.Fatalf("BubbleP returned an unexpected error: %v", err)
	}

	dew, err := raoult.DewP(raoult.SaturationPressureInput{
		Compositions: z, PSats: psat,
	})
	if err != nil {
		t.Fatalf("DewP returned an unexpected error: %v", err)
	}

	if !(dew.P < bubble.P) {
		t.Fatalf("the dew pressure %.4f does not lie below the bubble pressure %.4f", dew.P, bubble.P)
	}

	t.Run("at the bubble pressure the feed is all liquid", func(t *testing.T) {
		// Approached from inside the region: exactly at the boundary the
		// classification reports a single phase, which is correct.
		got, err := raoult.FlashPT(raoult.SaturationPressureInput{
			P: bubble.P * (1 - 1e-12), Compositions: z, PSats: psat, Options: tight(),
		})
		if err != nil {
			t.Fatalf("FlashPT returned an unexpected error: %v", err)
		}

		if math.Abs(got.V) > tol {
			t.Errorf("vapor fraction = %.9f at the bubble pressure; want 0", got.V)
		}

		// The liquid must still be the feed, having lost nothing to the
		// vapor.
		for i := range z {
			if math.Abs(got.X[i]-z[i]) > tol {
				t.Errorf("x%d = %.9f at the bubble pressure; want the feed %.4f", i+1, got.X[i], z[i])
			}
		}

		// And the vapor must be what a bubble-point calculation predicts.
		for i := range z {
			if math.Abs(got.Y[i]-bubble.Y[i]) > tol {
				t.Errorf("y%d = %.9f; the bubble-point calculation gives %.9f", i+1, got.Y[i], bubble.Y[i])
			}
		}
	})

	t.Run("at the dew pressure the feed is all vapor", func(t *testing.T) {
		got, err := raoult.FlashPT(raoult.SaturationPressureInput{
			P: dew.P * (1 + 1e-12), Compositions: z, PSats: psat, Options: tight(),
		})
		if err != nil {
			t.Fatalf("FlashPT returned an unexpected error: %v", err)
		}

		if math.Abs(got.V-1) > tol {
			t.Errorf("vapor fraction = %.9f at the dew pressure; want 1", got.V)
		}

		for i := range z {
			if math.Abs(got.Y[i]-z[i]) > tol {
				t.Errorf("y%d = %.9f at the dew pressure; want the feed %.4f", i+1, got.Y[i], z[i])
			}
		}

		for i := range z {
			if math.Abs(got.X[i]-dew.X[i]) > tol {
				t.Errorf("x%d = %.9f; the dew-point calculation gives %.9f", i+1, got.X[i], dew.X[i])
			}
		}
	})
}

// TestFlashMaterialBalance checks that the two phases account for the
// feed exactly, at every pressure across the two-phase region.
func TestFlashMaterialBalance(t *testing.T) {
	const tol = 1e-12

	z := []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	for _, pressure := range []float64{102, 105, 110, 120, 130} {
		got, err := raoult.FlashPT(raoult.SaturationPressureInput{
			P: pressure, Compositions: z, PSats: psat, Options: tight(),
		})
		if err != nil {
			t.Fatalf("FlashPT at %g kPa returned an unexpected error: %v", pressure, err)
		}

		for i := range z {
			balance := got.X[i]*got.L + got.Y[i]*got.V

			if math.Abs(balance-z[i]) > tol {
				t.Errorf(
					"at %g kPa, component %d: the phases hold %.15f; the feed is %.15f",
					pressure, i, balance, z[i],
				)
			}
		}
	}
}

// TestFlashEquilibriumRatios checks that the two phases stand in the
// ratio Raoult's law prescribes.
//
//	yi / xi = Pi_sat / P
func TestFlashEquilibriumRatios(t *testing.T) {
	const tol = 1e-9

	z := []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	const pressure = 110.0

	got, err := raoult.FlashPT(raoult.SaturationPressureInput{
		P: pressure, Compositions: z, PSats: psat, Options: tight(),
	})
	if err != nil {
		t.Fatalf("FlashPT returned an unexpected error: %v", err)
	}

	for i := range z {
		want := psat[i] / pressure

		if ratio := got.Y[i] / got.X[i]; math.Abs(ratio-want) > tol {
			t.Errorf("component %d: y/x = %.9f; want Psat/P = %.9f", i, ratio, want)
		}
	}
}

// TestFlashCompositionsAreNormalised checks that both phases are reported
// as compositions.
//
// They sum to one only when the vapor fraction is the exact root, so the
// tolerance here reflects the accuracy the solver was asked for.
func TestFlashCompositionsAreNormalised(t *testing.T) {
	const tol = 1e-9

	got, err := raoult.FlashPT(raoult.SaturationPressureInput{
		P:            110,
		Compositions: []float64{0.45, 0.35, 0.20},
		PSats:        []float64{195.75, 97.84, 50.32},
		Options:      tight(),
	})
	if err != nil {
		t.Fatalf("FlashPT returned an unexpected error: %v", err)
	}

	sum := func(v []float64) float64 {
		var total float64
		for _, e := range v {
			total += e
		}
		return total
	}

	if s := sum(got.X); math.Abs(s-1) > tol {
		t.Errorf("liquid composition sums to %.12f; want 1", s)
	}

	if s := sum(got.Y); math.Abs(s-1) > tol {
		t.Errorf("vapor composition sums to %.12f; want 1", s)
	}
}

// TestFlashOutsideTheTwoPhaseRegion checks that a feed which does not
// separate is described rather than reported as a failure.
func TestFlashOutsideTheTwoPhaseRegion(t *testing.T) {
	z := []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	testCases := []struct {
		name     string
		pressure float64
		want     vle.PhaseState
	}{
		{"compressed past the bubble point", 200, vle.SubcooledLiquid},
		{"expanded past the dew point", 60, vle.SuperheatedVapor},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			_, err := raoult.FlashPT(raoult.SaturationPressureInput{
				P: tc.pressure, Compositions: z, PSats: psat, Options: tight(),
			})

			if err == nil {
				t.Fatal("expected an error; got nil")
			}

			var single *vle.SinglePhaseError
			if !errors.As(err, &single) {
				t.Fatalf("error is not a *vle.SinglePhaseError: %v", err)
			}

			if single.State != tc.want {
				t.Errorf("state = %v; want %v", single.State, tc.want)
			}
		})
	}
}

// TestFlashBinaryAgainstTheExplicitSolution checks the general solver
// against the closed form a binary admits.
//
// With two components Raoult's law can be inverted directly:
//
//	x1 = (P - P2_sat) / (P1_sat - P2_sat)
//
// The general path solves the Rachford-Rice equation instead, so
// agreement is evidence that the iteration lands where the algebra says
// it should.
func TestFlashBinaryAgainstTheExplicitSolution(t *testing.T) {
	const tol = 1e-9

	psat := []float64{195.75, 50.32}
	z := []float64{0.5, 0.5}

	// The explicit inversion only describes a two-phase state between the
	// boundaries, which for a binary are available in closed form.
	bubble := z[0]*psat[0] + z[1]*psat[1]
	dew := 1 / (z[0]/psat[0] + z[1]/psat[1])

	for _, fraction := range []float64{0.1, 0.3, 0.5, 0.7, 0.9} {
		pressure := dew + fraction*(bubble-dew)

		got, err := raoult.FlashPT(raoult.SaturationPressureInput{
			P: pressure, Compositions: z, PSats: psat, Options: tight(),
		})
		if err != nil {
			t.Fatalf("FlashPT at %g kPa returned an unexpected error: %v", pressure, err)
		}

		want := (pressure - psat[1]) / (psat[0] - psat[1])

		if math.Abs(got.X[0]-want) > tol {
			t.Errorf("at %g kPa: x1 = %.12f; the explicit solution gives %.12f", pressure, got.X[0], want)
		}
	}
}

// TestFlashInvalidInput checks the guards on the feed and the conditions.
func TestFlashInvalidInput(t *testing.T) {
	psat := []float64{195.75, 97.84, 50.32}

	testCases := []struct {
		name  string
		input raoult.SaturationPressureInput
	}{
		{
			name: "composition does not sum to one",
			input: raoult.SaturationPressureInput{
				P: 110, Compositions: []float64{0.3, 0.3, 0.3}, PSats: psat,
			},
		},
		{
			name: "negative mole fraction",
			input: raoult.SaturationPressureInput{
				P: 110, Compositions: []float64{-0.2, 0.8, 0.4}, PSats: psat,
			},
		},
		{
			name: "more components than saturation pressures",
			input: raoult.SaturationPressureInput{
				P: 110, Compositions: []float64{0.25, 0.25, 0.25, 0.25}, PSats: psat,
			},
		},
		{
			name: "a single component cannot separate",
			input: raoult.SaturationPressureInput{
				P: 110, Compositions: []float64{1}, PSats: []float64{195.75},
			},
		},
		{
			name: "non-positive pressure",
			input: raoult.SaturationPressureInput{
				P: 0, Compositions: []float64{0.45, 0.35, 0.20}, PSats: psat,
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := raoult.FlashPT(tc.input); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}
