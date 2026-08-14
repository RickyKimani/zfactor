package internal

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/vle"
)

// tightOptions converges far below the accuracy any assertion here needs,
// so that a failure indicates a defect rather than a tolerance.
func tightOptions() vle.SolverOptions {
	return vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 200}
}

// example13_8 returns the feed and equilibrium ratios of Example 13.8 of
// Smith, Van Ness & Abbott: acetone(1)/acetonitrile(2)/nitromethane(3) at
// 80 °C and 110 kPa, where Raoult's law makes Ki = Pi_sat/P.
func example13_8() (z, k []float64) {
	z = []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	const pressure = 110.0

	k = make([]float64, len(psat))
	for i, p := range psat {
		k[i] = p / pressure
	}

	return z, k
}

// TestVaporFractionExample13_8 reproduces the worked flash calculation.
//
// The published answers are a vapor fraction of 0.7364, compositions
// y = {0.5087, 0.3389, 0.1524} and x = {0.2859, 0.3810, 0.3331}.
func TestVaporFractionExample13_8(t *testing.T) {
	z, k := example13_8()

	// The published equilibrium ratios, as a check on the state itself.
	for i, want := range []float64{1.7795, 0.8895, 0.4575} {
		if math.Abs(k[i]-want) > 5e-4 {
			t.Errorf("K%d = %.4f; want %.4f", i+1, k[i], want)
		}
	}

	v, err := VaporFraction(z, k, tightOptions())
	if err != nil {
		t.Fatalf("VaporFraction returned an unexpected error: %v", err)
	}

	if math.Abs(v-0.7364) > 5e-4 {
		t.Errorf("vapor fraction = %.6f; want 0.7364", v)
	}

	x, y, err := FlashCompositions(z, k, v)
	if err != nil {
		t.Fatalf("FlashCompositions returned an unexpected error: %v", err)
	}

	for i, want := range []float64{0.5087, 0.3389, 0.1524} {
		if math.Abs(y[i]-want) > 5e-4 {
			t.Errorf("y%d = %.4f; want %.4f", i+1, y[i], want)
		}
	}

	for i, want := range []float64{0.2859, 0.3810, 0.3331} {
		if math.Abs(x[i]-want) > 5e-4 {
			t.Errorf("x%d = %.4f; want %.4f", i+1, x[i], want)
		}
	}
}

// TestFlashCompositionsSatisfyTheMaterialBalance checks the identity the
// compositions are constructed from:
//
//	zi = xi(1 - V) + yi V.
//
// Substituting yi = zi Ki/(1 + V(Ki - 1)) and xi = yi/Ki reduces the
// right-hand side to zi for any vapor fraction whatever, so the balance
// holds exactly and independently of whether V solves the Rachford-Rice
// equation. The compositions summing to one is the part that depends on
// V being the root; that is checked separately below.
func TestFlashCompositionsSatisfyTheMaterialBalance(t *testing.T) {
	const tol = 1e-15

	z, k := example13_8()

	for _, v := range []float64{0, 0.1, 0.25, 0.5, 0.7364, 0.9, 1} {
		x, y, err := FlashCompositions(z, k, v)
		if err != nil {
			t.Fatalf("FlashCompositions returned an unexpected error: %v", err)
		}

		for i := range z {
			got := x[i]*(1-v) + y[i]*v

			if math.Abs(got-z[i]) > tol {
				t.Errorf(
					"at V = %g, component %d: x(1-V) + yV = %.15f; want the feed %.15f",
					v, i, got, z[i],
				)
			}
		}
	}
}

// TestFlashCompositionsSumToOneAtTheRoot checks the condition that
// singles out the solution.
//
// Both compositions sum to one only when the vapor fraction solves the
// Rachford-Rice equation, which is what the equation states. Away from
// the root they do not, and the departure is what the solver drives to
// zero.
func TestFlashCompositionsSumToOneAtTheRoot(t *testing.T) {
	z, k := example13_8()

	v, err := VaporFraction(z, k, tightOptions())
	if err != nil {
		t.Fatalf("VaporFraction returned an unexpected error: %v", err)
	}

	x, y, err := FlashCompositions(z, k, v)
	if err != nil {
		t.Fatalf("FlashCompositions returned an unexpected error: %v", err)
	}

	sum := func(v []float64) float64 {
		var total float64
		for _, e := range v {
			total += e
		}
		return total
	}

	if got := sum(x); math.Abs(got-1) > 1e-9 {
		t.Errorf("liquid composition sums to %.12f; want 1", got)
	}

	if got := sum(y); math.Abs(got-1) > 1e-9 {
		t.Errorf("vapor composition sums to %.12f; want 1", got)
	}

	// Away from the root they must not, or the root would carry no
	// information.
	offRoot := v / 2

	x, y, err = FlashCompositions(z, k, offRoot)
	if err != nil {
		t.Fatalf("FlashCompositions returned an unexpected error: %v", err)
	}

	if math.Abs(sum(x)-1) < 1e-6 && math.Abs(sum(y)-1) < 1e-6 {
		t.Errorf(
			"at V = %g, away from the root, both compositions still sum to one (%.9f, %.9f)",
			offRoot, sum(x), sum(y),
		)
	}
}

// TestRachfordRiceAvoidsTheTrivialRoots checks the reason for solving the
// difference of the two summation conditions rather than either alone.
//
// Summing the vapor compositions gives an identity at V = 1, and the
// liquid compositions at V = 0, so each single condition carries a root
// that is not a phase split. A solver applied to one of them can be drawn
// to that endpoint. Their difference vanishes at neither.
func TestRachfordRiceAvoidsTheTrivialRoots(t *testing.T) {
	z, k := example13_8()

	sumVapor := func(v float64) float64 {
		var total float64
		for i := range z {
			total += z[i] * k[i] / (1 + v*(k[i]-1))
		}
		return total - 1
	}

	sumLiquid := func(v float64) float64 {
		var total float64
		for i := range z {
			total += z[i] / (1 + v*(k[i]-1))
		}
		return total - 1
	}

	difference := func(v float64) float64 {
		var total float64
		for i := range z {
			total += z[i] * (k[i] - 1) / (1 + v*(k[i]-1))
		}
		return total
	}

	if got := sumVapor(1); math.Abs(got) > 1e-15 {
		t.Errorf("the vapor condition at V = 1 is %.3e; expected the trivial root", got)
	}

	if got := sumLiquid(0); math.Abs(got) > 1e-15 {
		t.Errorf("the liquid condition at V = 0 is %.3e; expected the trivial root", got)
	}

	if got := difference(0); math.Abs(got) < 1e-6 {
		t.Errorf("the difference vanishes at V = 0 (%.3e); it should not", got)
	}

	if got := difference(1); math.Abs(got) < 1e-6 {
		t.Errorf("the difference vanishes at V = 1 (%.3e); it should not", got)
	}
}

// TestRachfordRiceIsMonotonic checks the property that makes the root
// unique and the bracket reliable.
//
// The derivative of the difference form is a sum of squares with a
// leading minus sign, so it is negative wherever defined. The function
// therefore decreases across the whole interval, crossing zero once.
func TestRachfordRiceIsMonotonic(t *testing.T) {
	z, k := example13_8()

	difference := func(v float64) float64 {
		var total float64
		for i := range z {
			total += z[i] * (k[i] - 1) / (1 + v*(k[i]-1))
		}
		return total
	}

	previous := difference(0)

	for step := 1; step <= 100; step++ {
		v := float64(step) / 100
		got := difference(v)

		if got >= previous {
			t.Errorf("at V = %.2f the residual is %.9f; it should fall below the previous %.9f", v, got, previous)
		}

		previous = got
	}
}

// TestRachfordRicePolesLieOutsideTheInterval checks the structural fact
// that lets bisection run over the whole range without guarding against
// singularities.
//
// The poles sit at V = 1/(1 - Ki). For a ratio above one that is
// negative, and for one below it exceeds unity, so no pole falls within
// [0, 1] for any positive set of ratios.
func TestRachfordRicePolesLieOutsideTheInterval(t *testing.T) {
	for _, k := range [][]float64{
		{1.7795, 0.8895, 0.4575},
		{5, 0.01},
		{1.0001, 0.9999},
		{100, 50, 0.001, 0.5},
	} {
		for i, ki := range k {
			if ki == 1 {
				continue // no pole; the term vanishes
			}

			pole := 1 / (1 - ki)

			if pole >= 0 && pole <= 1 {
				t.Errorf("ratio %d of %v places a pole at V = %g, inside the interval", i, k, pole)
			}
		}
	}
}

// TestClassifyFeed checks the three states a feed can occupy, and that
// the boundaries agree with the bubble and dew conditions.
//
// A feed splits only when Σ zi Ki exceeds one, which places it below its
// bubble pressure, and Σ zi/Ki falls below one, which places it above its
// dew pressure.
func TestClassifyFeed(t *testing.T) {
	z := []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	ratios := func(pressure float64) []float64 {
		k := make([]float64, len(psat))
		for i, p := range psat {
			k[i] = p / pressure
		}
		return k
	}

	// Under Raoult's law the boundaries are the bubble and dew pressures
	// of the feed, which are available in closed form.
	var bubble, dewReciprocal float64
	for i := range z {
		bubble += z[i] * psat[i]
		dewReciprocal += z[i] / psat[i]
	}
	dew := 1 / dewReciprocal

	testCases := []struct {
		name     string
		pressure float64
		want     vle.PhaseState
	}{
		{"above the bubble pressure", bubble * 1.1, vle.SubcooledLiquid},
		{"just above the bubble pressure", bubble * 1.000001, vle.SubcooledLiquid},
		{"inside the two-phase region", (bubble + dew) / 2, vle.TwoPhase},
		{"just below the dew pressure", dew * 0.999999, vle.SuperheatedVapor},
		{"below the dew pressure", dew * 0.9, vle.SuperheatedVapor},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := ClassifyFeed(z, ratios(tc.pressure))
			if err != nil {
				t.Fatalf("ClassifyFeed returned an unexpected error: %v", err)
			}

			if got != tc.want {
				t.Errorf("at %.4f kPa the feed is %v; want %v", tc.pressure, got, tc.want)
			}
		})
	}
}

// TestVaporFractionSinglePhase checks that a feed which does not split is
// reported as such, with the state recoverable from the error.
func TestVaporFractionSinglePhase(t *testing.T) {
	z := []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	ratios := func(pressure float64) []float64 {
		k := make([]float64, len(psat))
		for i, p := range psat {
			k[i] = p / pressure
		}
		return k
	}

	testCases := []struct {
		name     string
		pressure float64
		want     vle.PhaseState
	}{
		{"subcooled liquid", 200, vle.SubcooledLiquid},
		{"superheated vapor", 60, vle.SuperheatedVapor},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := VaporFraction(z, ratios(tc.pressure), tightOptions())

			if err == nil {
				t.Fatalf("expected an error; got a vapor fraction of %g", got)
			}

			var single *vle.SinglePhaseError
			if !errors.As(err, &single) {
				t.Fatalf("error is not a *vle.SinglePhaseError: %v", err)
			}

			if single.State != tc.want {
				t.Errorf("state = %v; want %v", single.State, tc.want)
			}

			if got != 0 {
				t.Errorf("vapor fraction = %g; want 0 alongside the error", got)
			}
		})
	}
}

// TestVaporFractionDecreasesWithPressure checks the direction of the
// split: compressing a two-phase mixture condenses it.
func TestVaporFractionDecreasesWithPressure(t *testing.T) {
	z := []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	previous := math.Inf(1)

	for _, pressure := range []float64{105, 110, 115, 120, 125, 130} {
		k := make([]float64, len(psat))
		for i, p := range psat {
			k[i] = p / pressure
		}

		v, err := VaporFraction(z, k, tightOptions())
		if err != nil {
			t.Fatalf("VaporFraction at %g kPa returned an unexpected error: %v", pressure, err)
		}

		if v >= previous {
			t.Errorf("at %g kPa the vapor fraction is %.6f; it should fall below the previous %.6f", pressure, v, previous)
		}

		previous = v
	}
}

// TestVaporFractionInvalidInput checks the guards on the feed and the
// equilibrium ratios.
//
// A non-positive ratio is rejected because it would place a pole of the
// residual inside the search interval, and because it describes a
// component present in one phase and wholly absent from the other.
func TestVaporFractionInvalidInput(t *testing.T) {
	testCases := []struct {
		name string
		z, k []float64
	}{
		{"no components", nil, nil},
		{"mismatched lengths", []float64{0.5, 0.5}, []float64{2}},
		{"composition does not sum to one", []float64{0.3, 0.3}, []float64{2, 0.5}},
		{"negative mole fraction", []float64{-0.2, 1.2}, []float64{2, 0.5}},
		{"zero equilibrium ratio", []float64{0.5, 0.5}, []float64{2, 0}},
		{"negative equilibrium ratio", []float64{0.5, 0.5}, []float64{2, -0.5}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := VaporFraction(tc.z, tc.k, tightOptions()); err == nil {
				t.Error("VaporFraction: expected an error; got nil")
			}

			if _, _, err := FlashCompositions(tc.z, tc.k, 0.5); err == nil {
				t.Error("FlashCompositions: expected an error; got nil")
			}
		})
	}
}

// TestFlashCompositionsRejectsImpossibleFraction checks that a vapor
// fraction outside the physical range is refused.
func TestFlashCompositionsRejectsImpossibleFraction(t *testing.T) {
	z, k := example13_8()

	for _, v := range []float64{-0.1, 1.1, 2} {
		if _, _, err := FlashCompositions(z, k, v); err == nil {
			t.Errorf("at V = %g: expected an error; got nil", v)
		}
	}
}

// TestPhaseStateString checks that each state renders a name, since the
// single-phase error embeds one in its message.
func TestPhaseStateString(t *testing.T) {
	testCases := map[vle.PhaseState]string{
		vle.TwoPhase:         "two-phase",
		vle.SubcooledLiquid:  "subcooled liquid",
		vle.SuperheatedVapor: "superheated vapor",
		vle.PhaseState(99):   "unknown",
	}

	for state, want := range testCases {
		if got := state.String(); got != want {
			t.Errorf("state %d renders as %q; want %q", int(state), got, want)
		}
	}
}
