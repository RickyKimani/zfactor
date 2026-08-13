package internal

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/vle"
)

// errSentinel stands in for a failure raised by the residual itself, so
// that error propagation can be distinguished from a solver's own
// failure modes.
var errSentinel = errors.New("residual failed")

// noError adapts a plain function to the signature the solvers accept.
func noError(f func(float64) float64) func(float64) (float64, error) {
	return func(x float64) (float64, error) {
		return f(x), nil
	}
}

// TestBisectionFindsKnownRoots checks the solver against functions whose
// roots are known in closed form.
//
// Bisection converges for any continuous function that changes sign
// across the bracket, so the cases include a polynomial, a transcendental
// function and one with a root at an endpoint of the search.
func TestBisectionFindsKnownRoots(t *testing.T) {
	const tol = 1e-10

	opts := vle.SolverOptions{Tolerance: tol, MaxIterations: 200}

	testCases := []struct {
		name   string
		f      func(float64) float64
		lo, hi float64
		want   float64
	}{
		{"square root of two", func(x float64) float64 { return x*x - 2 }, 0, 2, math.Sqrt2},
		{"cube root of ten", func(x float64) float64 { return x*x*x - 10 }, 0, 5, math.Cbrt(10)},
		{"cosine fixed point", func(x float64) float64 { return math.Cos(x) - x }, 0, 1, 0.7390851332151607},
		{"natural logarithm", func(x float64) float64 { return math.Log(x) }, 0.1, 10, 1},
		{"descending sign change", func(x float64) float64 { return 2 - x }, 0, 5, 2},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := Bisection(noError(tc.f), tc.lo, tc.hi, opts)
			if err != nil {
				t.Fatalf("Bisection returned an unexpected error: %v", err)
			}

			if math.Abs(got-tc.want) > tol {
				t.Errorf("root = %.12f; want %.12f", got, tc.want)
			}

			// The iteration must never leave the interval it was given.
			if got < tc.lo || got > tc.hi {
				t.Errorf("root %.12f lies outside the bracket [%g, %g]", got, tc.lo, tc.hi)
			}
		})
	}
}

// TestBisectionStaysWithinBracket checks the property that motivates
// using bisection in place of the secant method.
//
// The azeotrope search relies on it: a mole fraction must remain inside
// (0, 1), and an open-ended method can step outside a domain where the
// residual is undefined. Here the residual reports an error anywhere
// beyond the bracket, so any excursion is caught rather than merely
// producing a poor estimate.
func TestBisectionStaysWithinBracket(t *testing.T) {
	const lo, hi = 0.0, 1.0

	guarded := func(x float64) (float64, error) {
		if x < lo || x > hi {
			t.Errorf("solver evaluated the residual at %g, outside [%g, %g]", x, lo, hi)
			return 0, errSentinel
		}
		return x - 0.25, nil
	}

	got, err := Bisection(guarded, lo, hi, vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 200})
	if err != nil {
		t.Fatalf("Bisection returned an unexpected error: %v", err)
	}

	if math.Abs(got-0.25) > 1e-12 {
		t.Errorf("root = %.15f; want 0.25", got)
	}
}

// TestBisectionRootAtEndpoint checks that a residual vanishing exactly at
// a bound is returned directly rather than bisected toward.
func TestBisectionRootAtEndpoint(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 100}

	t.Run("lower bound", func(t *testing.T) {
		got, err := Bisection(noError(func(x float64) float64 { return x }), 0, 3, opts)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if got != 0 {
			t.Errorf("root = %g; want exactly 0", got)
		}
	})

	t.Run("upper bound", func(t *testing.T) {
		got, err := Bisection(noError(func(x float64) float64 { return x - 3 }), 0, 3, opts)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}
		if got != 3 {
			t.Errorf("root = %g; want exactly 3", got)
		}
	})
}

// TestBisectionRejectsBadBracket checks that an interval which does not
// enclose a root is refused rather than returning a number from it.
//
// Bisection has no way to locate a root it does not bracket, so silently
// returning the midpoint of a same-sign interval would be worse than
// failing.
func TestBisectionRejectsBadBracket(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-9, MaxIterations: 100}

	testCases := []struct {
		name   string
		f      func(float64) float64
		lo, hi float64
	}{
		{"no sign change, both positive", func(x float64) float64 { return x*x + 1 }, 0, 1},
		{"no sign change, both negative", func(x float64) float64 { return -x*x - 1 }, 0, 1},
		{"two roots inside, same sign at the ends", func(x float64) float64 { return (x - 1) * (x - 2) }, 0, 3},
		{"bounds reversed", func(x float64) float64 { return x }, 1, -1},
		{"empty interval", func(x float64) float64 { return x }, 1, 1},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := Bisection(noError(tc.f), tc.lo, tc.hi, opts); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestBisectionRespectsIterationLimit checks that an exhausted iteration
// budget is reported rather than a half-converged value returned.
func TestBisectionRespectsIterationLimit(t *testing.T) {
	// One iteration cannot narrow a wide bracket to this tolerance.
	opts := vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 1}

	if _, err := Bisection(noError(func(x float64) float64 { return x - 0.3 }), 0, 1, opts); err == nil {
		t.Error("expected a convergence error; got nil")
	}
}

// TestBisectionPropagatesResidualError checks that a failure inside the
// residual reaches the caller unchanged, rather than being reported as a
// convergence failure.
func TestBisectionPropagatesResidualError(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-9, MaxIterations: 100}

	failAt := func(threshold float64) func(float64) (float64, error) {
		return func(x float64) (float64, error) {
			if x > threshold {
				return 0, errSentinel
			}
			return x - 0.5, nil
		}
	}

	// Fails when first evaluating the upper bound.
	if _, err := Bisection(failAt(0.9), 0, 1, opts); !errors.Is(err, errSentinel) {
		t.Errorf("error = %v; want the sentinel from the residual", err)
	}

	// Fails partway through, at a midpoint.
	if _, err := Bisection(failAt(0.6), 0, 1, opts); !errors.Is(err, errSentinel) {
		t.Errorf("error = %v; want the sentinel from the residual", err)
	}
}

// TestSecantFindsKnownRoots checks the secant method against roots known
// in closed form.
func TestSecantFindsKnownRoots(t *testing.T) {
	const tol = 1e-10

	opts := vle.SolverOptions{Tolerance: tol, MaxIterations: 100}

	testCases := []struct {
		name   string
		f      func(float64) float64
		x0, x1 float64
		want   float64
	}{
		{"square root of two", func(x float64) float64 { return x*x - 2 }, 1, 2, math.Sqrt2},
		{"cosine fixed point", func(x float64) float64 { return math.Cos(x) - x }, 0, 1, 0.7390851332151607},
		{"linear, exact in one step", func(x float64) float64 { return 3*x - 6 }, 0, 1, 2},
		{"exponential", func(x float64) float64 { return math.Exp(x) - 5 }, 1, 2, math.Log(5)},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := Secant(noError(tc.f), tc.x0, tc.x1, opts)
			if err != nil {
				t.Fatalf("Secant returned an unexpected error: %v", err)
			}

			if math.Abs(got-tc.want) > 1e-8 {
				t.Errorf("root = %.12f; want %.12f", got, tc.want)
			}
		})
	}
}

// TestSecantRejectsFlatResidual checks the guard against a vanishing
// slope.
//
// The secant step divides by the difference between two residual values.
// A constant residual makes that difference zero, which would otherwise
// produce an infinity or a NaN and carry it into the result.
func TestSecantRejectsFlatResidual(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-9, MaxIterations: 100}

	if _, err := Secant(noError(func(float64) float64 { return 1 }), 0, 1, opts); err == nil {
		t.Error("expected an error for a constant residual; got nil")
	}
}

// TestSecantRespectsIterationLimit checks that an exhausted iteration
// budget is reported rather than a partly converged value returned.
func TestSecantRespectsIterationLimit(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 2}

	// The root is at sqrt(2); two steps from these distant guesses close
	// only a fraction of the gap.
	if _, err := Secant(noError(func(x float64) float64 { return x*x - 2 }), 100, 101, opts); err == nil {
		t.Error("expected a convergence error; got nil")
	}
}

// TestSecantConvergesOnStepSizeAlone documents a limitation of the
// convergence criterion rather than asserting a desirable outcome.
//
// Secant stops once successive iterates differ by less than the
// tolerance, without checking that the residual has actually reached
// zero. Where the function is very flat the step collapses first, and
// the method reports success at a point that need not be a root.
//
// Here exp(x) has no root at all, yet the iteration settles at x = -50
// and returns without error, the residual being about 2e-22 there. For
// the VLE residuals this is harmless — they are steep where they are
// used — but a caller supplying its own function should know that a
// returned value is not by itself evidence of a root.
func TestSecantConvergesOnStepSizeAlone(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-15, MaxIterations: 100}

	got, err := Secant(noError(math.Exp), -50, 50, opts)
	if err != nil {
		t.Skipf("behaviour has changed; Secant now reports %v", err)
	}

	if residual := math.Exp(got); residual > 1e-9 {
		t.Errorf(
			"Secant returned %g where the residual is %g; the criterion is expected to stop only where the residual is negligible",
			got, residual,
		)
	}
}

// TestSecantPropagatesResidualError checks that a failure inside the
// residual reaches the caller, from both the initial evaluations and the
// iteration.
func TestSecantPropagatesResidualError(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-9, MaxIterations: 100}

	alwaysFails := func(float64) (float64, error) { return 0, errSentinel }

	if _, err := Secant(alwaysFails, 0, 1, opts); !errors.Is(err, errSentinel) {
		t.Errorf("error = %v; want the sentinel from the residual", err)
	}

	calls := 0
	failsLater := func(x float64) (float64, error) {
		calls++
		if calls > 2 {
			return 0, errSentinel
		}
		return x*x - 2, nil
	}

	if _, err := Secant(failsLater, 0, 1, opts); !errors.Is(err, errSentinel) {
		t.Errorf("error = %v; want the sentinel from the residual", err)
	}
}

// TestFixedPointConverges checks successive substitution on maps that
// contract, where the iteration is guaranteed to converge to the unique
// fixed point.
func TestFixedPointConverges(t *testing.T) {
	const tol = 1e-12

	opts := vle.SolverOptions{Tolerance: tol, MaxIterations: 500}

	t.Run("cosine map", func(t *testing.T) {
		// x = cos(x) has the Dottie number as its only fixed point.
		got, err := FixedPoint(
			[]float64{1},
			func(x []float64) ([]float64, error) {
				return []float64{math.Cos(x[0])}, nil
			},
			opts,
		)
		if err != nil {
			t.Fatalf("FixedPoint returned an unexpected error: %v", err)
		}

		if math.Abs(got[0]-0.7390851332151607) > 1e-9 {
			t.Errorf("fixed point = %.12f; want 0.739085133215", got[0])
		}
	})

	t.Run("coupled pair", func(t *testing.T) {
		// Each component halves toward a distinct target, so the
		// iteration exercises the vector convergence criterion rather
		// than a scalar one.
		got, err := FixedPoint(
			[]float64{0, 0},
			func(x []float64) ([]float64, error) {
				return []float64{(x[0] + 4) / 2, (x[1] + 10) / 2}, nil
			},
			opts,
		)
		if err != nil {
			t.Fatalf("FixedPoint returned an unexpected error: %v", err)
		}

		if math.Abs(got[0]-4) > 1e-9 || math.Abs(got[1]-10) > 1e-9 {
			t.Errorf("fixed point = %v; want [4 10]", got)
		}
	})

	t.Run("already at the fixed point", func(t *testing.T) {
		got, err := FixedPoint(
			[]float64{2},
			func(x []float64) ([]float64, error) { return []float64{x[0]}, nil },
			opts,
		)
		if err != nil {
			t.Fatalf("FixedPoint returned an unexpected error: %v", err)
		}

		if got[0] != 2 {
			t.Errorf("fixed point = %g; want 2", got[0])
		}
	})
}

// TestFixedPointDoesNotModifyInput checks that the caller's slice is left
// untouched, since the solver copies it before iterating.
func TestFixedPointDoesNotModifyInput(t *testing.T) {
	initial := []float64{1, 2}
	want := []float64{1, 2}

	_, err := FixedPoint(
		initial,
		func(x []float64) ([]float64, error) {
			return []float64{(x[0] + 4) / 2, (x[1] + 10) / 2}, nil
		},
		vle.SolverOptions{Tolerance: 1e-12, MaxIterations: 500},
	)
	if err != nil {
		t.Fatalf("FixedPoint returned an unexpected error: %v", err)
	}

	for i := range initial {
		if initial[i] != want[i] {
			t.Errorf("input was modified: got %v, want %v", initial, want)
			break
		}
	}
}

// TestFixedPointDiverges checks that a map which does not contract is
// reported as a failure rather than returning whatever the last iterate
// happened to be.
func TestFixedPointDiverges(t *testing.T) {
	opts := vle.SolverOptions{Tolerance: 1e-9, MaxIterations: 50}

	_, err := FixedPoint(
		[]float64{1},
		func(x []float64) ([]float64, error) {
			return []float64{2 * x[0]}, nil
		},
		opts,
	)

	if err == nil {
		t.Error("expected a convergence error for a divergent map; got nil")
	}
}

// TestFixedPointPropagatesUpdateError checks that a failure inside the
// update reaches the caller.
func TestFixedPointPropagatesUpdateError(t *testing.T) {
	_, err := FixedPoint(
		[]float64{1},
		func([]float64) ([]float64, error) { return nil, errSentinel },
		vle.SolverOptions{Tolerance: 1e-9, MaxIterations: 50},
	)

	if !errors.Is(err, errSentinel) {
		t.Errorf("error = %v; want the sentinel from the update", err)
	}
}

// TestSolverOptionDefaults checks that a zero-valued SolverOptions is
// usable, as its documentation promises, by solving with one.
func TestSolverOptionDefaults(t *testing.T) {
	var opts vle.SolverOptions

	if opts.Tol() <= 0 {
		t.Errorf("default tolerance = %g; want a positive value", opts.Tol())
	}

	if opts.MaxIter() <= 0 {
		t.Errorf("default iteration limit = %d; want a positive value", opts.MaxIter())
	}

	got, err := Bisection(noError(func(x float64) float64 { return x*x - 2 }), 0, 2, opts)
	if err != nil {
		t.Fatalf("Bisection with default options returned an unexpected error: %v", err)
	}

	if math.Abs(got-math.Sqrt2) > opts.Tol() {
		t.Errorf("root = %.12f; want %.12f", got, math.Sqrt2)
	}
}
