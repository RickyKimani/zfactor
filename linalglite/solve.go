package linalglite

import (
	"errors"
	"fmt"
	"math"
)

// Solve returns the x satisfying A·x = b.
//
// It is the short form, for a system solved once. An iteration that
// solves repeatedly against the same matrix should keep a factorisation
// from Factorize instead, and one that changes the matrix each step
// should reuse it through Refactorize.
//
// Systems of order three and below are solved by their closed forms,
// which cost less than a factorisation and allocate nothing beyond the
// result. Those orders cover a good deal of this field: a Newton step on
// a ternary equilibrium, a two-parameter fit, the extents of a pair of
// simultaneous reactions.
func Solve(a *Dense, b []float64) ([]float64, error) {
	if !a.IsSquare() {
		return nil, fmt.Errorf(
			"%w: solving needs a square matrix, got %dx%d",
			ErrShape, a.rows, a.cols,
		)
	}

	if len(b) != a.rows {
		return nil, fmt.Errorf(
			"%w: an order-%d system needs %d right-hand side values, got %d",
			ErrShape, a.rows, a.rows, len(b),
		)
	}

	if x, done, err := solveSmall(a, b); done {
		return x, err
	}

	lu, err := Factorize(a)
	if err != nil {
		return nil, err
	}

	return lu.Solve(b)
}

// solveSmall handles the orders that have a closed form, reporting
// whether it dealt with the system.
//
// The forms are Cramer's rule, which for these sizes is both faster than
// elimination and free of the pivot search. It is only worth using while
// the number of terms stays small: the determinant of an order-n matrix
// expands into n! of them, so this stops at three.
func solveSmall(a *Dense, b []float64) (x []float64, done bool, err error) {
	switch a.rows {
	case 1:
		if a.data[0] == 0 {
			return nil, true, fmt.Errorf("%w: the single element is zero", ErrSingular)
		}

		return []float64{b[0] / a.data[0]}, true, nil

	case 2:
		m := a.data

		det := m[0]*m[3] - m[1]*m[2]
		if det == 0 {
			return nil, true, fmt.Errorf("%w: the determinant is zero", ErrSingular)
		}

		inverse := 1 / det

		return []float64{
			(b[0]*m[3] - b[1]*m[1]) * inverse,
			(b[1]*m[0] - b[0]*m[2]) * inverse,
		}, true, nil

	case 3:
		m := a.data

		// The cofactors of the first row, reused for the determinant and
		// for each of the three numerators.
		c00 := m[4]*m[8] - m[5]*m[7]
		c01 := m[5]*m[6] - m[3]*m[8]
		c02 := m[3]*m[7] - m[4]*m[6]

		det := m[0]*c00 + m[1]*c01 + m[2]*c02
		if det == 0 {
			return nil, true, fmt.Errorf("%w: the determinant is zero", ErrSingular)
		}

		inverse := 1 / det

		// Cramer's rule with the substituted columns expanded, so that
		// each unknown costs a handful of multiplications.
		return []float64{
			(b[0]*c00 + b[1]*(m[2]*m[7]-m[1]*m[8]) + b[2]*(m[1]*m[5]-m[2]*m[4])) * inverse,
			(b[0]*c01 + b[1]*(m[0]*m[8]-m[2]*m[6]) + b[2]*(m[2]*m[3]-m[0]*m[5])) * inverse,
			(b[0]*c02 + b[1]*(m[1]*m[6]-m[0]*m[7]) + b[2]*(m[0]*m[4]-m[1]*m[3])) * inverse,
		}, true, nil
	}

	return nil, false, nil
}

// Det returns the determinant of a square matrix.
//
// Orders three and below use their closed forms; larger matrices are
// factorised, the determinant being the product of the diagonal of the
// upper factor with a sign from the row interchanges.
//
// A singular matrix has a determinant of zero, which is a value rather
// than an error, so this returns one where Solve would report ErrSingular.
func Det(a *Dense) (float64, error) {
	if !a.IsSquare() {
		return 0, fmt.Errorf(
			"%w: a determinant needs a square matrix, got %dx%d",
			ErrShape, a.rows, a.cols,
		)
	}

	m := a.data

	switch a.rows {
	case 1:
		return m[0], nil
	case 2:
		return m[0]*m[3] - m[1]*m[2], nil
	case 3:
		return m[0]*(m[4]*m[8]-m[5]*m[7]) -
			m[1]*(m[3]*m[8]-m[5]*m[6]) +
			m[2]*(m[3]*m[7]-m[4]*m[6]), nil
	}

	lu, err := Factorize(a)
	if err != nil {
		// A singular matrix has a determinant, and it is zero.
		if errors.Is(err, ErrSingular) {
			return 0, nil
		}

		return 0, err
	}

	return lu.Det(), nil
}

// Residual returns the largest absolute element of A·x - b, the direct
// measure of how well x solves the system.
//
// It is offered because the useful check on a solution is not how it was
// obtained but whether it satisfies the equations. Scaling it by the
// magnitudes involved gives a relative figure comparable across problems.
func Residual(a *Dense, x, b []float64) (float64, error) {
	product, err := a.MulVec(x)
	if err != nil {
		return 0, err
	}

	if len(b) != len(product) {
		return 0, fmt.Errorf(
			"%w: the product has %d elements and the right-hand side %d",
			ErrShape, len(product), len(b),
		)
	}

	var worst float64
	for i := range product {
		if d := math.Abs(product[i] - b[i]); d > worst {
			worst = d
		}
	}

	return worst, nil
}
