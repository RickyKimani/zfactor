package internal

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor/vle"
)

// Secant solves f(x) = 0 using the secant method.
//
// The method requires two initial guesses x0 and x1.
//
// Convergence is achieved when:
//
//	|x(k+1) - x(k)| < tolerance
//
// The criterion is on the step rather than on the residual. Since the
// step is the residual scaled by the secant slope, a small step does
// imply a residual small against the local variation of the function,
// which is the appropriate comparison for a solver that knows nothing
// about the scale of what it is given. Where a function is very flat,
// though, the step can collapse while the residual remains meaningful,
// so a returned value is not on its own proof of a root.
//
// The solver returns an error if the residual takes the same value at
// both points, if a step leaves the finite range, or if it fails to
// converge within the maximum iteration limit.
func Secant(
	f func(float64) (float64, error),
	x0, x1 float64,
	opts vle.SolverOptions,
) (float64, error) {

	tol := opts.Tol()
	maxIter := opts.MaxIter()

	f0, err := f(x0)
	if err != nil {
		return 0, err
	}

	f1, err := f(x1)
	if err != nil {
		return 0, err
	}

	for range maxIter {

		denom := f1 - f0

		// The secant step divides by the change in the residual. Only an
		// exact zero makes that impossible; comparing against a fixed
		// magnitude instead would reject well-posed problems whose
		// residuals are simply small, since the threshold would then sit
		// above the whole scale of the function.
		if denom == 0 {
			return 0, errors.New(
				"secant method failed: the residual takes the same value at both points",
			)
		}

		x2 := x1 - f1*(x1-x0)/denom

		if math.IsNaN(x2) || math.IsInf(x2, 0) {
			return 0, errors.New(
				"secant method failed: the step left the finite range",
			)
		}

		if math.Abs(x2-x1) < tol {
			return x2, nil
		}

		x0 = x1
		f0 = f1

		x1 = x2

		f1, err = f(x1)
		if err != nil {
			return 0, err
		}
	}

	return 0, errors.New(
		"secant method failed to converge",
	)
}

// Bisection solves f(x) = 0 on the bracketing interval [lo, hi].
//
// The method requires the residual to change sign across the interval,
// which guarantees that a root is enclosed and that the iteration
// cannot leave the bracket. This makes it the appropriate choice when
// the search must stay inside a physically meaningful domain, such as a
// mole fraction confined to (0, 1), where the open-ended secant method
// may wander outside.
//
// Convergence is achieved when the half-width of the bracket falls
// below the specified tolerance.
//
// The solver returns an error if the interval is malformed or if the
// residual does not change sign across it.
func Bisection(
	f func(float64) (float64, error),
	lo, hi float64,
	opts vle.SolverOptions,
) (float64, error) {

	if lo >= hi {
		return 0, errors.New(
			"bisection method failed: lower bound must be below upper bound",
		)
	}

	tol := opts.Tol()
	maxIter := opts.MaxIter()

	flo, err := f(lo)
	if err != nil {
		return 0, err
	}
	if flo == 0 {
		return lo, nil
	}

	fhi, err := f(hi)
	if err != nil {
		return 0, err
	}
	if fhi == 0 {
		return hi, nil
	}

	if math.Signbit(flo) == math.Signbit(fhi) {
		return 0, errors.New(
			"bisection method failed: residual does not change sign across the interval",
		)
	}

	for range maxIter {
		mid := (lo + hi) / 2

		if (hi-lo)/2 < tol {
			return mid, nil
		}

		fmid, err := f(mid)
		if err != nil {
			return 0, err
		}
		if fmid == 0 {
			return mid, nil
		}

		if math.Signbit(fmid) == math.Signbit(flo) {
			lo, flo = mid, fmid
		} else {
			hi = mid
		}
	}

	return 0, errors.New("bisection method failed to converge")
}

// FixedPoint solves x = g(x) by successive substitution.
//
// Convergence is achieved when the maximum absolute change between
// successive iterates is less than the specified tolerance.
func FixedPoint(
	initial []float64,
	update func([]float64) ([]float64, error),
	opts vle.SolverOptions,
) ([]float64, error) {

	tol := opts.Tol()
	maxIter := opts.MaxIter()

	x := make([]float64, len(initial))
	copy(x, initial)

	for range maxIter {
		next, err := update(x)
		if err != nil {
			return nil, err
		}

		maxDiff := 0.0
		for i := range x {
			diff := math.Abs(next[i] - x[i])
			if diff > maxDiff {
				maxDiff = diff
			}
		}

		if maxDiff < tol {
			return next, nil
		}

		x = next
	}

	return nil, errors.New("fixed-point iteration failed to converge")
}
