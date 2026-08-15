package linalglite

import (
	"fmt"
	"math"
)

// LU is the LU factorisation of a square matrix with partial pivoting.
//
// The factors are held in one buffer, the unit-diagonal lower triangle
// below the diagonal and the upper triangle on and above it, together
// with the row interchanges that produced them. Solving against the
// factorisation costs a fraction of computing it, so a matrix used
// against several right-hand sides should be factorised once.
//
// An LU is not safe for concurrent use: SolveInto reads the factors and
// writes the destination, and Refactorize overwrites the factors
// entirely.
type LU struct {
	n     int
	lu    []float64 // the factors, row-major, n by n
	pivot []int     // pivot[k] is the row swapped into position k
	swaps int       // number of interchanges, for the sign of the determinant

	// scratch holds the permuted right-hand side during a solve, so that
	// SolveInto need not allocate and need not modify its input.
	scratch []float64
}

// Factorize computes the LU factorisation of a square matrix.
//
// The matrix is not modified; its elements are copied into storage the
// factorisation owns. Where a loop factorises a succession of matrices of
// the same size, Refactorize reuses that storage.
//
// It returns ErrSingular if the matrix has no unique inverse, and ErrShape
// if it is not square.
func Factorize(a *Dense) (*LU, error) {
	if !a.IsSquare() {
		return nil, fmt.Errorf(
			"%w: factorisation needs a square matrix, got %dx%d",
			ErrShape, a.rows, a.cols,
		)
	}

	n := a.rows

	lu := &LU{
		n:       n,
		lu:      make([]float64, n*n),
		pivot:   make([]int, n),
		scratch: make([]float64, n),
	}

	if err := lu.Refactorize(a); err != nil {
		return nil, err
	}

	return lu, nil
}

// Refactorize replaces the factorisation with that of a, reusing the
// storage already held.
//
// It is the allocation-free path for iterations that change the matrix
// each step, such as a Newton method on a Jacobian. The matrix must have
// the same order as the one originally factorised.
func (lu *LU) Refactorize(a *Dense) error {
	if !a.IsSquare() {
		return fmt.Errorf(
			"%w: factorisation needs a square matrix, got %dx%d",
			ErrShape, a.rows, a.cols,
		)
	}

	if a.rows != lu.n {
		return fmt.Errorf(
			"%w: this factorisation holds an order-%d matrix, got order %d",
			ErrShape, lu.n, a.rows,
		)
	}

	copy(lu.lu, a.data)

	return lu.decompose()
}

// decompose performs the factorisation in place on lu.lu.
//
// The algorithm is right-looking Doolittle with partial pivoting: for
// each column in turn, the largest remaining entry is brought to the
// diagonal, the entries below it are divided by it to form that column of
// the lower factor, and the trailing submatrix is updated by the outer
// product of the column and row just fixed.
//
// The pivot search is what keeps the result accurate. Eliminating with a
// small pivot multiplies the rest of the matrix by its reciprocal, and
// the precision lost in that step is not recovered by any later one.
func (lu *LU) decompose() error {
	n := lu.n
	data := lu.lu

	lu.swaps = 0

	for k := range n {
		// Find the largest entry in the remaining part of column k.
		pivot := k
		largest := math.Abs(data[k*n+k])

		for i := k + 1; i < n; i++ {
			if v := math.Abs(data[i*n+k]); v > largest {
				largest, pivot = v, i
			}
		}

		if largest == 0 {
			// The column is entirely zero below the diagonal, so no
			// interchange can supply a pivot: the rows are dependent.
			return fmt.Errorf("%w: no pivot in column %d", ErrSingular, k)
		}

		lu.pivot[k] = pivot

		if pivot != k {
			// Swap the two rows whole. Physically moving them keeps every
			// later access contiguous, which is worth more than the
			// indirection an index permutation would save.
			rowK := data[k*n : k*n+n]
			rowP := data[pivot*n : pivot*n+n]

			for j := range rowK {
				rowK[j], rowP[j] = rowP[j], rowK[j]
			}

			lu.swaps++
		}

		rowK := data[k*n : k*n+n]
		diagonal := rowK[k]

		// Reciprocal once, then multiply: division is several times the
		// cost of multiplication, and this runs n times per column.
		inverse := 1 / diagonal

		// The part of row k that updates the trailing submatrix.
		above := rowK[k+1 : n]

		for i := k + 1; i < n; i++ {
			rowI := data[i*n : i*n+n]

			factor := rowI[k] * inverse
			rowI[k] = factor

			if factor == 0 {
				// Nothing to subtract. Stoichiometric and Jacobian
				// matrices in this field are often sparse enough for this
				// to skip real work.
				continue
			}

			// Re-slicing to the same length as above tells the compiler
			// the two walks have equal extent, so the indexed access
			// needs no bounds check.
			trailing := rowI[k+1 : n]
			trailing = trailing[:len(above)]

			for j, v := range above {
				trailing[j] -= factor * v
			}
		}
	}

	return nil
}

// Order returns the order of the factorised matrix.
func (lu *LU) Order() int { return lu.n }

// Det returns the determinant of the factorised matrix.
//
// It is the product of the diagonal of the upper factor, negated once for
// each row interchange. Computing it from the factorisation costs n
// multiplications, against the n! terms of the direct expansion.
//
// The product can overflow or underflow for a large matrix even when the
// determinant is representable; where only the sign or the magnitude
// matters, LogAbsDet avoids that.
func (lu *LU) Det() float64 {
	det := 1.0

	for i := range lu.n {
		det *= lu.lu[i*lu.n+i]
	}

	if lu.swaps%2 == 1 {
		return -det
	}

	return det
}

// LogAbsDet returns the natural logarithm of the absolute value of the
// determinant, together with its sign as -1, 0 or +1.
//
// Summing logarithms rather than multiplying diagonal entries keeps the
// result representable for matrices whose determinant is not.
func (lu *LU) LogAbsDet() (log float64, sign float64) {
	sign = 1
	if lu.swaps%2 == 1 {
		sign = -1
	}

	for i := range lu.n {
		v := lu.lu[i*lu.n+i]

		if v == 0 {
			return math.Inf(-1), 0
		}

		if v < 0 {
			sign = -sign
			v = -v
		}

		log += math.Log(v)
	}

	return log, sign
}

// Solve returns the x satisfying A·x = b for the factorised A.
func (lu *LU) Solve(b []float64) ([]float64, error) {
	x := make([]float64, len(b))

	if err := lu.SolveInto(x, b); err != nil {
		return nil, err
	}

	return x, nil
}

// SolveInto writes the solution of A·x = b into dst, allocating nothing.
//
// dst must hold one element per equation. It may be the same slice as b,
// which is convenient for an iteration that overwrites its residual with
// the step.
func (lu *LU) SolveInto(dst, b []float64) error {
	n := lu.n

	if len(b) != n {
		return fmt.Errorf(
			"%w: an order-%d system needs %d right-hand side values, got %d",
			ErrShape, n, n, len(b),
		)
	}

	if len(dst) != n {
		return fmt.Errorf(
			"%w: an order-%d system produces %d values, but the destination holds %d",
			ErrShape, n, n, len(dst),
		)
	}

	// Apply the row interchanges to a copy, leaving b untouched and
	// allowing dst and b to be the same slice.
	x := lu.scratch
	copy(x, b)

	for k := range n {
		if p := lu.pivot[k]; p != k {
			x[k], x[p] = x[p], x[k]
		}
	}

	data := lu.lu

	// Forward substitution through the unit lower triangle. Its diagonal
	// is implicitly one, so no division is needed here.
	for i := 1; i < n; i++ {
		row := data[i*n : i*n+i]

		var sum float64
		for j, v := range row {
			sum += v * x[j]
		}

		x[i] -= sum
	}

	// Back substitution through the upper triangle.
	for i := n - 1; i >= 0; i-- {
		row := data[i*n+i+1 : i*n+n]
		tail := x[i+1 : n]
		tail = tail[:len(row)]

		var sum float64
		for j, v := range row {
			sum += v * tail[j]
		}

		x[i] = (x[i] - sum) / data[i*n+i]
	}

	copy(dst, x)

	return nil
}

// SolveMatrixInto writes the solution of A·X = B into dst, one column of
// B at a time.
//
// dst and b must have as many rows as the system has equations and the
// same number of columns as each other. They may be the same matrix.
func (lu *LU) SolveMatrixInto(dst, b *Dense) error {
	if b.rows != lu.n || dst.rows != lu.n {
		return fmt.Errorf(
			"%w: an order-%d system needs %d rows, got %d and %d",
			ErrShape, lu.n, lu.n, b.rows, dst.rows,
		)
	}

	if dst.cols != b.cols {
		return fmt.Errorf(
			"%w: the destination has %d columns and the right-hand side %d",
			ErrShape, dst.cols, b.cols,
		)
	}

	column := make([]float64, lu.n)

	for j := range b.cols {
		for i := range lu.n {
			column[i] = b.data[i*b.cols+j]
		}

		if err := lu.SolveInto(column, column); err != nil {
			return err
		}

		for i := range lu.n {
			dst.data[i*dst.cols+j] = column[i]
		}
	}

	return nil
}

// Inverse returns the inverse of the factorised matrix.
//
// Solving A·x = b directly is both faster and more accurate than forming
// the inverse and multiplying by it, so this is for the cases that genuinely
// want the matrix itself.
func (lu *LU) Inverse() (*Dense, error) {
	inverse := Identity(lu.n)

	if err := lu.SolveMatrixInto(inverse, inverse); err != nil {
		return nil, err
	}

	return inverse, nil
}
