package raoult

import (
	"errors"
	"math"
)

const (
	defaultRootTolerance = 1e-6
	defaultMaxIterations = 10
)

type SolverOptions struct {
	Tolerance     float64
	MaxIterations int
}

func (s SolverOptions) tolerance() float64 {
	if s.Tolerance <= 0 {
		return defaultRootTolerance
	}
	return s.Tolerance
}

func (s SolverOptions) maxIterations() int {
	if s.MaxIterations <= 0 {
		return defaultMaxIterations
	}
	return s.MaxIterations
}

// secant solves f(x) = 0 using the secant method.
//
// The method requires two initial guesses x0 and x1.
//
// Convergence is achieved when:
//
//	|x(k+1) - x(k)| < tolerance
//
// The solver returns an error if the method encounters a near-zero
// slope or fails to converge within the maximum iteration limit.
func secant(
	f func(float64) (float64, error),
	x0, x1 float64,
	opts SolverOptions,
) (float64, error) {

	tol := opts.tolerance()
	maxIter := opts.maxIterations()

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
