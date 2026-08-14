package linalglite

import (
	"errors"
	"fmt"
)

// ErrShape reports a matrix or vector whose size does not fit the
// operation asked of it.
var ErrShape = errors.New("linalglite: shape mismatch")

// ErrSingular reports a matrix with no unique inverse.
//
// Factorisation detects this when a column has no non-zero entry left to
// pivot on, meaning the rows are linearly dependent. A system with no
// unique solution is reported rather than answered, since dividing by the
// vanishing pivot would return infinities that look like numbers.
var ErrSingular = errors.New("linalglite: matrix is singular")

// Dense is a matrix stored by rows in one contiguous slice.
//
// The zero value is not usable; construct with New, NewFrom or Identity.
type Dense struct {
	rows, cols int
	data       []float64
}

// New returns a rows-by-cols matrix of zeros.
//
// It panics if either dimension is not positive, since a matrix with no
// elements is a programming error rather than a runtime condition.
func New(rows, cols int) *Dense {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("linalglite: New(%d, %d): dimensions must be positive", rows, cols))
	}

	return &Dense{
		rows: rows,
		cols: cols,
		data: make([]float64, rows*cols),
	}
}

// NewFrom returns a rows-by-cols matrix holding the given elements, read
// in row-major order.
//
// The slice is copied, so the caller may reuse or modify it afterwards.
// It panics if the length does not match the dimensions.
func NewFrom(rows, cols int, data []float64) *Dense {
	if rows <= 0 || cols <= 0 {
		panic(fmt.Sprintf("linalglite: NewFrom(%d, %d): dimensions must be positive", rows, cols))
	}

	if len(data) != rows*cols {
		panic(fmt.Sprintf(
			"linalglite: NewFrom(%d, %d): got %d elements, want %d",
			rows, cols, len(data), rows*cols,
		))
	}

	m := &Dense{
		rows: rows,
		cols: cols,
		data: make([]float64, len(data)),
	}
	copy(m.data, data)

	return m
}

// Identity returns the n-by-n matrix with ones on its diagonal.
func Identity(n int) *Dense {
	m := New(n, n)

	for i := range n {
		m.data[i*n+i] = 1
	}

	return m
}

// Rows returns the number of rows.
func (m *Dense) Rows() int { return m.rows }

// Cols returns the number of columns.
func (m *Dense) Cols() int { return m.cols }

// IsSquare reports whether the matrix has as many rows as columns.
func (m *Dense) IsSquare() bool { return m.rows == m.cols }

// At returns the element in row i and column j, counting from zero.
//
// It panics if either index is out of range.
func (m *Dense) At(i, j int) float64 {
	if uint(i) >= uint(m.rows) || uint(j) >= uint(m.cols) {
		panic(fmt.Sprintf("linalglite: At(%d, %d) on a %dx%d matrix", i, j, m.rows, m.cols))
	}

	return m.data[i*m.cols+j]
}

// Set writes v into row i and column j, counting from zero.
//
// It panics if either index is out of range.
func (m *Dense) Set(i, j int, v float64) {
	if uint(i) >= uint(m.rows) || uint(j) >= uint(m.cols) {
		panic(fmt.Sprintf("linalglite: Set(%d, %d) on a %dx%d matrix", i, j, m.rows, m.cols))
	}

	m.data[i*m.cols+j] = v
}

// Row returns the elements of row i as a slice into the matrix.
//
// Writing through it changes the matrix. It is offered so that a caller
// filling a Jacobian row by row need not go through Set for every
// element, which is the difference between one bounds check and n.
func (m *Dense) Row(i int) []float64 {
	if uint(i) >= uint(m.rows) {
		panic(fmt.Sprintf("linalglite: Row(%d) on a %dx%d matrix", i, m.rows, m.cols))
	}

	return m.data[i*m.cols : (i+1)*m.cols : (i+1)*m.cols]
}

// RawRowMajor returns the underlying slice, in row-major order.
//
// Writing through it changes the matrix. It exists for bulk filling and
// for interoperating with code that works in flat slices.
func (m *Dense) RawRowMajor() []float64 { return m.data }

// Clone returns a copy holding the same elements.
func (m *Dense) Clone() *Dense {
	c := &Dense{
		rows: m.rows,
		cols: m.cols,
		data: make([]float64, len(m.data)),
	}
	copy(c.data, m.data)

	return c
}

// Zero sets every element to zero, keeping the storage.
//
// Filling a fresh Jacobian each iteration is the usual reason to want
// this; it avoids the allocation a new matrix would cost.
func (m *Dense) Zero() {
	clear(m.data)
}

// MulVec returns the product of the matrix with the vector x.
//
// It is provided mainly so that a caller can check a solution: the
// residual A·x - b is the direct measure of whether a solve succeeded.
func (m *Dense) MulVec(x []float64) ([]float64, error) {
	dst := make([]float64, m.rows)

	if err := m.MulVecInto(dst, x); err != nil {
		return nil, err
	}

	return dst, nil
}

// MulVecInto writes the product of the matrix with x into dst,
// allocating nothing.
//
// dst must have one element per row and must not alias x.
func (m *Dense) MulVecInto(dst, x []float64) error {
	if len(x) != m.cols {
		return fmt.Errorf(
			"%w: cannot multiply a %dx%d matrix by a vector of length %d",
			ErrShape, m.rows, m.cols, len(x),
		)
	}

	if len(dst) != m.rows {
		return fmt.Errorf(
			"%w: a %dx%d matrix produces %d elements, but the destination holds %d",
			ErrShape, m.rows, m.cols, m.rows, len(dst),
		)
	}

	for i := range m.rows {
		row := m.data[i*m.cols : i*m.cols+m.cols]

		// Both slices are length m.cols, which lets the compiler drop the
		// bounds check on the indexed one.
		xs := x[:len(row)]

		var sum float64
		for j, v := range row {
			sum += v * xs[j]
		}

		dst[i] = sum
	}

	return nil
}
