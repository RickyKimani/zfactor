package modified_raoult

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/vle"
	"github.com/rickykimani/zfactor/vle/raoult"
)

// tightFlash converges far below the accuracy any assertion here needs,
// so a failure indicates a defect rather than a tolerance.
func tightFlash() vle.SolverOptions {
	return vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 500}
}

// flashInput builds a feed for the methanol/methyl acetate system of
// Example 13.1 at 318.15 K.
func flashInput(pressure float64, z []float64, activity ActivityModel) MixtureInput {
	return MixtureInput{
		T:            exampleTemperature,
		P:            pressure,
		Compositions: z,
		Antoine:      exampleModels(),
		Activity:     activity,
		Options:      tightFlash(),
	}
}

// twoPhaseWindow returns the pressures between which a feed separates,
// found from the bubble- and dew-point calculations.
func twoPhaseWindow(t *testing.T, z []float64, activity ActivityModel) (dew, bubble float64) {
	t.Helper()

	b, err := BubbleP(flashInput(0, z, activity))
	if err != nil {
		t.Fatalf("BubbleP returned an unexpected error: %v", err)
	}

	d, err := DewP(flashInput(0, z, activity))
	if err != nil {
		t.Fatalf("DewP returned an unexpected error: %v", err)
	}

	if !(d.P < b.P) {
		t.Fatalf("the dew pressure %.4f does not lie below the bubble pressure %.4f", d.P, b.P)
	}

	return d.P, b.P
}

// TestFlashPTIdealActivityMatchesRaoult checks that the non-ideal flash
// collapses onto the ideal one when the solution is ideal.
//
// Margules parameters of zero give unit activity coefficients at every
// composition, at which point the equilibrium ratios reduce to
// Pi_sat/P and the iteration has nothing left to do. The two packages
// reach that state by different routes — one solves the Rachford-Rice
// equation once, the other wraps it in a loop over the activity
// coefficients and precedes it with two boundary calculations — so
// agreement is evidence that the loop converges to the right place.
func TestFlashPTIdealActivityMatchesRaoult(t *testing.T) {
	const tol = 1e-9

	ideal := Margules{A12: 0, A21: 0}
	z := []float64{0.4, 0.6}

	dew, bubble := twoPhaseWindow(t, z, ideal)

	for _, fraction := range []float64{0.1, 0.25, 0.5, 0.75, 0.9} {
		pressure := dew + fraction*(bubble-dew)

		modified, err := FlashPT(flashInput(pressure, z, ideal))
		if err != nil {
			t.Fatalf("FlashPT at %.4f kPa returned an unexpected error: %v", pressure, err)
		}

		simple, err := raoult.FlashPT(raoult.SaturationPressureInput{
			P:            pressure,
			Compositions: z,
			PSats:        []float64{methanolPSat(t), methylAcetatePSat(t)},
			Options:      tightFlash(),
		})
		if err != nil {
			t.Fatalf("raoult.FlashPT at %.4f kPa returned an unexpected error: %v", pressure, err)
		}

		if math.Abs(modified.V-simple.V) > tol {
			t.Errorf(
				"at %.4f kPa: vapor fraction %.12f against Raoult's %.12f",
				pressure, modified.V, simple.V,
			)
		}

		for i := range z {
			if math.Abs(modified.X[i]-simple.X[i]) > tol {
				t.Errorf("at %.4f kPa: x%d = %.12f against Raoult's %.12f", pressure, i+1, modified.X[i], simple.X[i])
			}

			if math.Abs(modified.Y[i]-simple.Y[i]) > tol {
				t.Errorf("at %.4f kPa: y%d = %.12f against Raoult's %.12f", pressure, i+1, modified.Y[i], simple.Y[i])
			}
		}
	}
}

// methanolPSat and methylAcetatePSat evaluate the example's correlations
// at its temperature, for comparison against the ideal-law package.
func methanolPSat(t *testing.T) float64 {
	t.Helper()

	p, err := methanolAntoine().Pressure(exampleTemperature)
	if err != nil {
		t.Fatalf("Pressure returned an unexpected error: %v", err)
	}

	return p
}

func methylAcetatePSat(t *testing.T) float64 {
	t.Helper()

	p, err := methylAcetateAntoine().Pressure(exampleTemperature)
	if err != nil {
		t.Fatalf("Pressure returned an unexpected error: %v", err)
	}

	return p
}

// TestFlashPTMaterialBalance checks that the two phases account for the
// feed exactly.
//
// The compositions are constructed from the material balance, so it holds
// however well or badly the outer iteration has converged. It is
// therefore a check on the construction rather than on the convergence.
func TestFlashPTMaterialBalance(t *testing.T) {
	const tol = 1e-12

	activity := exampleActivity()
	z := []float64{0.5, 0.5}

	dew, bubble := twoPhaseWindow(t, z, activity)

	for _, fraction := range []float64{0.05, 0.25, 0.5, 0.75, 0.95} {
		pressure := dew + fraction*(bubble-dew)

		got, err := FlashPT(flashInput(pressure, z, activity))
		if err != nil {
			t.Fatalf("FlashPT at %.4f kPa returned an unexpected error: %v", pressure, err)
		}

		for i := range z {
			balance := got.X[i]*got.L + got.Y[i]*got.V

			if math.Abs(balance-z[i]) > tol {
				t.Errorf(
					"at %.4f kPa, component %d: the phases hold %.15f; the feed is %.15f",
					pressure, i, balance, z[i],
				)
			}
		}

		if math.Abs(got.V+got.L-1) > tol {
			t.Errorf("at %.4f kPa the phases hold %.15f mol between them; want 1", pressure, got.V+got.L)
		}
	}
}

// TestFlashPTSatisfiesTheEquilibriumRelation checks the condition the
// converged state must obey.
//
// The activity coefficients are evaluated at the liquid composition the
// flash reports, and the partial pressures must then agree across the
// phase boundary:
//
//	yi P = xi γi Pi_sat.
//
// This is the modified Raoult's law itself, so it is the definition of a
// converged answer rather than a value read from anywhere.
func TestFlashPTSatisfiesTheEquilibriumRelation(t *testing.T) {
	const relTol = 1e-6

	activity := exampleActivity()
	z := []float64{0.45, 0.55}

	dew, bubble := twoPhaseWindow(t, z, activity)

	psat, err := MixtureInput{T: exampleTemperature, Antoine: exampleModels()}.PSat()
	if err != nil {
		t.Fatalf("PSat returned an unexpected error: %v", err)
	}

	for _, fraction := range []float64{0.2, 0.5, 0.8} {
		pressure := dew + fraction*(bubble-dew)

		got, err := FlashPT(flashInput(pressure, z, activity))
		if err != nil {
			t.Fatalf("FlashPT at %.4f kPa returned an unexpected error: %v", pressure, err)
		}

		gamma, err := activityCoefficients(activity, exampleTemperature, got.X)
		if err != nil {
			t.Fatalf("activityCoefficients returned an unexpected error: %v", err)
		}

		for i := range z {
			vapor := got.Y[i] * pressure
			liquid := got.X[i] * gamma[i] * psat[i]

			scale := math.Max(math.Abs(vapor), 1)

			if rel := math.Abs(vapor-liquid) / scale; rel > relTol {
				t.Errorf(
					"at %.4f kPa, component %d: yP = %.9f but xγPsat = %.9f (%.2e apart, relative)",
					pressure, i, vapor, liquid, rel,
				)
			}
		}
	}
}

// TestFlashPTApproachesTheBoundaries checks that the vapor fraction runs
// from nothing at the bubble pressure to everything at the dew pressure.
//
// Those are the two edges of the region, where the flash must agree with
// the calculations that located them.
func TestFlashPTApproachesTheBoundaries(t *testing.T) {
	const tol = 1e-3

	activity := exampleActivity()
	z := []float64{0.5, 0.5}

	dew, bubble := twoPhaseWindow(t, z, activity)

	nearBubble, err := FlashPT(flashInput(bubble*(1-1e-9), z, activity))
	if err != nil {
		t.Fatalf("FlashPT near the bubble pressure returned an unexpected error: %v", err)
	}

	if nearBubble.V > tol {
		t.Errorf("just inside the bubble pressure the vapor fraction is %.6f; want approximately 0", nearBubble.V)
	}

	nearDew, err := FlashPT(flashInput(dew*(1+1e-9), z, activity))
	if err != nil {
		t.Fatalf("FlashPT near the dew pressure returned an unexpected error: %v", err)
	}

	if nearDew.V < 1-tol {
		t.Errorf("just inside the dew pressure the vapor fraction is %.6f; want approximately 1", nearDew.V)
	}
}

// TestFlashPTVaporFractionFallsWithPressure checks the direction of the
// split: compressing a two-phase mixture condenses it.
func TestFlashPTVaporFractionFallsWithPressure(t *testing.T) {
	activity := exampleActivity()
	z := []float64{0.5, 0.5}

	dew, bubble := twoPhaseWindow(t, z, activity)

	previous := math.Inf(1)

	for _, fraction := range []float64{0.05, 0.2, 0.4, 0.6, 0.8, 0.95} {
		pressure := dew + fraction*(bubble-dew)

		got, err := FlashPT(flashInput(pressure, z, activity))
		if err != nil {
			t.Fatalf("FlashPT at %.4f kPa returned an unexpected error: %v", pressure, err)
		}

		if got.V >= previous {
			t.Errorf(
				"at %.4f kPa the vapor fraction is %.6f; it should fall below the previous %.6f",
				pressure, got.V, previous,
			)
		}

		previous = got.V
	}
}

// TestFlashPTOutsideTheTwoPhaseRegion checks that a feed which does not
// separate is described rather than reported as a failure.
//
// The bounds come from the bubble- and dew-point calculations rather than
// from the equilibrium ratios at the feed composition. Those ratios
// depend on the liquid composition, which near the dew point differs
// greatly from the feed, so reading the bounds from them would misplace
// the region.
func TestFlashPTOutsideTheTwoPhaseRegion(t *testing.T) {
	activity := exampleActivity()
	z := []float64{0.5, 0.5}

	dew, bubble := twoPhaseWindow(t, z, activity)

	testCases := []struct {
		name     string
		pressure float64
		want     vle.PhaseState
	}{
		{"compressed past the bubble point", bubble * 1.05, vle.SubcooledLiquid},
		{"expanded past the dew point", dew * 0.95, vle.SuperheatedVapor},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			_, err := FlashPT(flashInput(tc.pressure, z, activity))

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

// TestFlashPTWholeWindowSolves checks that every pressure between the
// boundaries yields a split.
//
// A gap would mean the region the flash accepts disagrees with the one
// the boundary calculations describe, which is how an earlier version
// failed: classifying from the activity coefficients at the feed
// composition rejected pressures near the dew point that do separate.
func TestFlashPTWholeWindowSolves(t *testing.T) {
	activity := exampleActivity()

	for _, z := range [][]float64{
		{0.25, 0.75},
		{0.5, 0.5},
		{0.8, 0.2},
	} {
		dew, bubble := twoPhaseWindow(t, z, activity)

		for step := 1; step < 20; step++ {
			pressure := dew + (bubble-dew)*float64(step)/20

			got, err := FlashPT(flashInput(pressure, z, activity))
			if err != nil {
				t.Errorf("feed %v at %.4f kPa, inside [%.4f, %.4f]: %v", z, pressure, dew, bubble, err)
				continue
			}

			if got.V < 0 || got.V > 1 {
				t.Errorf("feed %v at %.4f kPa: vapor fraction %.6f lies outside [0, 1]", z, pressure, got.V)
			}
		}
	}
}

// TestFlashPTInvalidInput checks the guards on the feed and the
// conditions.
func TestFlashPTInvalidInput(t *testing.T) {
	activity := exampleActivity()

	testCases := []struct {
		name  string
		input MixtureInput
	}{
		{
			name:  "composition does not sum to one",
			input: flashInput(70, []float64{0.3, 0.3}, activity),
		},
		{
			name:  "negative mole fraction",
			input: flashInput(70, []float64{-0.2, 1.2}, activity),
		},
		{
			name:  "non-positive pressure",
			input: flashInput(0, []float64{0.5, 0.5}, activity),
		},
		{
			name:  "no activity model",
			input: flashInput(70, []float64{0.5, 0.5}, nil),
		},
		{
			name: "more components than correlations",
			input: MixtureInput{
				T: exampleTemperature, P: 70,
				Compositions: []float64{0.3, 0.3, 0.4},
				Antoine:      exampleModels(),
				Activity:     activity,
				Options:      tightFlash(),
			},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := FlashPT(tc.input); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestFlashPTUsesTheSuppliedCorrelations checks that the Antoine models
// reach the calculation, by comparing against saturation pressures
// evaluated directly.
func TestFlashPTUsesTheSuppliedCorrelations(t *testing.T) {
	var _ []antoine.Model = exampleModels()

	activity := exampleActivity()
	z := []float64{0.5, 0.5}

	dew, bubble := twoPhaseWindow(t, z, activity)
	pressure := (dew + bubble) / 2

	got, err := FlashPT(flashInput(pressure, z, activity))
	if err != nil {
		t.Fatalf("FlashPT returned an unexpected error: %v", err)
	}

	gamma, err := activityCoefficients(activity, exampleTemperature, got.X)
	if err != nil {
		t.Fatalf("activityCoefficients returned an unexpected error: %v", err)
	}

	// Recover each saturation pressure from the converged state and
	// compare with the correlation.
	for i, want := range []float64{methanolPSat(t), methylAcetatePSat(t)} {
		recovered := got.Y[i] * pressure / (got.X[i] * gamma[i])

		if rel := math.Abs(recovered-want) / want; rel > 1e-6 {
			t.Errorf(
				"component %d: the converged state implies a saturation pressure of %.6f; the correlation gives %.6f",
				i, recovered, want,
			)
		}
	}
}
