package cubic

import (
	"math"
	"testing"
)

// soaveClone is an equation of state that behaves exactly like SRK but does
// not advertise its alpha derivative, so lnAlphaDeriv has to differentiate it
// numerically.
//
// Its purpose is to make the fallback comparable against a known answer: the
// closed form for the same alpha is available on the real SRK, so the two
// paths can be checked against each other over a range of states.
type soaveClone struct{}

func (soaveClone) Alpha(tr, w float64) float64 {
	b := 1 - math.Sqrt(tr)
	c := 1 + srkM(w)*b
	return c * c
}

func (soaveClone) Params() *Params {
	return &Params{Sigma: 1, Epsilon: 0, Omega: 0.08664, Psi: 0.42748}
}

// An equation that does not provide its derivative must still get a correct
// one. This differentiates a copy of SRK numerically and compares it against
// the closed form on SRK itself, which is the same function of the same alpha.
func TestLnAlphaDerivFallsBackToNumericalDifferentiation(t *testing.T) {
	var (
		exact    = &SRK{}
		fallback = soaveClone{}
	)

	// Confirm the premise, or the test would be comparing two closed forms
	// and the fallback would never run.
	if _, ok := interface{}(fallback).(LnAlphaDeriver); ok {
		t.Fatal("soaveClone advertises a derivative, so the fallback is not being exercised")
	}
	if _, ok := interface{}(exact).(LnAlphaDeriver); !ok {
		t.Fatal("SRK no longer provides a derivative")
	}

	for _, w := range []float64{0, 0.1, 0.3, 0.5} {
		for _, tr := range []float64{0.5, 0.7, 0.9, 1.0, 1.3, 2.0, 5.0} {
			want := lnAlphaDeriv(exact, tr, w)
			got := lnAlphaDeriv(fallback, tr, w)

			// The central difference carries the error here, not the
			// closed form, so the tolerance is set by the step size.
			const tolerance = 1e-7

			if diff := math.Abs(got - want); diff > tolerance*math.Max(1, math.Abs(want)) {
				t.Errorf("w=%g Tr=%g: numerical %.12g, closed form %.12g, differ by %.3g",
					w, tr, got, want, diff)
			}
		}
	}
}

// The alpha of the Soave form vanishes at sqrt(Tr) = 1 + 1/m, where ln alpha
// and so its derivative are undefined. The guard against it can only fire when
// the zero is exactly representable, which m = 1 arranges: sqrt(4) is 2
// exactly, so u is exactly zero there.
//
// For a real substance this lies far outside any useful range — for ethane
// under SRK it is above Tr = 4, some 1300 K — so this exercises the guard
// rather than a state anyone will reach.
func TestSoaveLnAlphaDerivIsNaNWhereAlphaVanishesExactly(t *testing.T) {
	const (
		m  = 1.0
		tr = 4.0
	)

	// Confirm the zero really is exact, or the guard would not be reached
	// and this would be testing the ordinary path.
	if u := 1 + m*(1-math.Sqrt(tr)); u != 0 {
		t.Fatalf("u is %g rather than exactly 0, so the guard is not being exercised", u)
	}

	if got := soaveLnAlphaDeriv(m, tr); !math.IsNaN(got) {
		t.Errorf("got %v where alpha vanishes, want NaN", got)
	}
}

// Approaching that zero, the derivative grows without bound: alpha tends to
// zero while its slope does not, so d ln alpha/d ln Tr diverges. This is the
// honest behaviour near the singularity, and it is worth pinning because the
// alternative — quietly clamping to a finite number — would hide a state the
// Soave form does not describe.
func TestSoaveLnAlphaDerivDivergesApproachingTheAlphaZero(t *testing.T) {
	const m = 1.0

	// Tr = 4 is the zero; these approach it from below.
	var previous float64

	for i, tr := range []float64{3.9, 3.99, 3.999, 3.9999} {
		got := math.Abs(soaveLnAlphaDeriv(m, tr))

		if math.IsNaN(got) || math.IsInf(got, 0) {
			t.Fatalf("Tr=%g: got %v short of the zero", tr, got)
		}

		if i > 0 && got <= previous {
			t.Errorf("Tr=%g: magnitude %.6g did not exceed %.6g from the previous step",
				tr, got, previous)
		}

		previous = got
	}
}

// A non-positive reduced temperature is not a state, and the square root in
// the Soave form would make it a NaN in any case.
func TestSoaveLnAlphaDerivRejectsNonPositiveTr(t *testing.T) {
	for _, tr := range []float64{0, -1} {
		if got := soaveLnAlphaDeriv(0.5, tr); !math.IsNaN(got) {
			t.Errorf("Tr=%g: got %v, want NaN", tr, got)
		}
	}
}

// The degenerate branch of the integral, taken when sigma equals epsilon, is
// the limit of the general one rather than a separate formula. Approaching
// that limit with an equation whose sigma and epsilon nearly coincide should
// reproduce the van der Waals value, which is what pins the branch as correct
// rather than merely finite.
func TestQTimesIDegenerateBranchIsTheLimitOfTheGeneral(t *testing.T) {
	const (
		Z = 0.9
		A = 0.1
		B = 0.05
	)

	vdw := &EOSCfg{Type: &VdW{}, T: 400, Tc: 305.32, Acentric: 0}

	degenerate, err := qTimesI(vdw, Z, A, B)
	if err != nil {
		t.Fatalf("van der Waals: %v", err)
	}

	// A sequence of equations whose sigma approaches epsilon from above.
	for _, gap := range []float64{1e-2, 1e-4, 1e-6} {
		cfg := &EOSCfg{Type: nearlyDegenerate{gap: gap}, T: 400, Tc: 305.32}

		general, err := qTimesI(cfg, Z, A, B)
		if err != nil {
			t.Fatalf("gap %g: %v", gap, err)
		}

		// The general form differs from the limit at first order in the
		// gap, so the agreement should tighten as the gap closes.
		if diff := math.Abs(general - degenerate); diff > gap {
			t.Errorf("gap %g: general form gives %.12g, the limit gives %.12g, differ by %.3g",
				gap, general, degenerate, diff)
		}
	}
}

// nearlyDegenerate is van der Waals with sigma pushed just off epsilon, so
// the general branch of the integral is taken for a case whose answer is
// known from the degenerate one.
type nearlyDegenerate struct{ gap float64 }

func (n nearlyDegenerate) Alpha(tr, w float64) float64 { return 1 }

func (n nearlyDegenerate) Params() *Params {
	return &Params{Sigma: n.gap, Epsilon: 0, Omega: 1.0 / 8.0, Psi: 27.0 / 64.0}
}
