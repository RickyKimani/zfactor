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
// The solver returns an error if the method encounters a near-zero
// slope or fails to converge within the maximum iteration limit.
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

		if math.Abs(denom) < 1e-14 {
			return 0, errors.New(
				"secant method failed: slope too close to zero",
			)
		}

		x2 := x1 - f1*(x1-x0)/denom

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
