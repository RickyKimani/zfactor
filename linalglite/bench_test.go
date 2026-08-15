package linalglite_test

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/rickykimani/zfactor/linalglite"
)

// benchMatrix returns a well-conditioned matrix of the given order.
//
// The diagonal is loaded so that no benchmark run stumbles into a
// singular matrix, which would measure the error path instead of the
// work.
func benchMatrix(n int) *linalglite.Dense {
	rng := rand.New(rand.NewSource(int64(n)))

	m := linalglite.New(n, n)

	for i := range n {
		row := m.Row(i)
		for j := range row {
			row[j] = rng.Float64()*2 - 1
		}
		row[i] += float64(n)
	}

	return m
}

func benchVector(n int) []float64 {
	rng := rand.New(rand.NewSource(int64(n) + 1))

	v := make([]float64, n)
	for i := range v {
		v[i] = rng.Float64()*2 - 1
	}

	return v
}

// BenchmarkSolve measures the short form, which factorises and discards.
//
// The orders below three take the closed forms; the rest factorise. The
// step between three and four is where that switch shows.
func BenchmarkSolve(b *testing.B) {
	for _, n := range []int{1, 2, 3, 4, 8, 16, 32, 64, 128} {
		a := benchMatrix(n)
		rhs := benchVector(n)

		b.Run(order(n), func(b *testing.B) {
			b.ReportAllocs()

			for b.Loop() {
				if _, err := linalglite.Solve(a, rhs); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// BenchmarkFactorize measures the factorisation alone, which is the cubic
// part of the work.
func BenchmarkFactorize(b *testing.B) {
	for _, n := range []int{4, 8, 16, 32, 64, 128} {
		a := benchMatrix(n)

		b.Run(order(n), func(b *testing.B) {
			b.ReportAllocs()

			for b.Loop() {
				if _, err := linalglite.Factorize(a); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// BenchmarkRefactorize measures replacing a factorisation in storage
// already held, which is the path a Newton iteration takes.
//
// Against BenchmarkFactorize it shows what the allocation costs.
func BenchmarkRefactorize(b *testing.B) {
	for _, n := range []int{4, 8, 16, 32, 64, 128} {
		a := benchMatrix(n)

		lu, err := linalglite.Factorize(a)
		if err != nil {
			b.Fatal(err)
		}

		b.Run(order(n), func(b *testing.B) {
			b.ReportAllocs()

			for b.Loop() {
				if err := lu.Refactorize(a); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// BenchmarkSolveInto measures a solve against a factorisation already
// computed, which is the quadratic part.
//
// This is what an iteration repeats when the matrix does not change, and
// against BenchmarkSolve it shows what keeping the factorisation saves.
func BenchmarkSolveInto(b *testing.B) {
	for _, n := range []int{4, 8, 16, 32, 64, 128} {
		a := benchMatrix(n)
		rhs := benchVector(n)

		lu, err := linalglite.Factorize(a)
		if err != nil {
			b.Fatal(err)
		}

		dst := make([]float64, n)

		b.Run(order(n), func(b *testing.B) {
			b.ReportAllocs()

			for b.Loop() {
				if err := lu.SolveInto(dst, rhs); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// BenchmarkNewtonStep measures the shape of work a simulation actually
// does: refill a Jacobian, refactorise it, solve for the step.
//
// It is the figure to watch, since it is the inner loop of every
// equilibrium and reaction calculation that will sit on this package.
func BenchmarkNewtonStep(b *testing.B) {
	for _, n := range []int{2, 3, 5, 10, 20, 50} {
		a := benchMatrix(n)
		residual := benchVector(n)

		lu, err := linalglite.Factorize(a)
		if err != nil {
			b.Fatal(err)
		}

		step := make([]float64, n)
		rng := rand.New(rand.NewSource(7))

		b.Run(order(n), func(b *testing.B) {
			b.ReportAllocs()

			for b.Loop() {
				// Stand in for evaluating the Jacobian afresh.
				for i := range n {
					row := a.Row(i)
					for j := range row {
						row[j] = rng.Float64()*2 - 1
					}
					row[i] += float64(n)
				}

				if err := lu.Refactorize(a); err != nil {
					b.Fatal(err)
				}

				if err := lu.SolveInto(step, residual); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// BenchmarkDet measures the determinant, whose small orders are closed
// forms and whose larger ones factorise.
func BenchmarkDet(b *testing.B) {
	for _, n := range []int{2, 3, 4, 16, 64} {
		a := benchMatrix(n)

		b.Run(order(n), func(b *testing.B) {
			b.ReportAllocs()

			for b.Loop() {
				if _, err := linalglite.Det(a); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// BenchmarkMulVecInto measures forming a residual, which an iteration
// does once per step alongside the solve.
func BenchmarkMulVecInto(b *testing.B) {
	for _, n := range []int{4, 16, 64, 128} {
		a := benchMatrix(n)
		x := benchVector(n)
		dst := make([]float64, n)

		b.Run(order(n), func(b *testing.B) {
			b.ReportAllocs()

			for b.Loop() {
				if err := a.MulVecInto(dst, x); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}

// order names a benchmark case by the order of its matrix, zero-padded so
// the output sorts readably.
func order(n int) string {
	return fmt.Sprintf("n=%03d", n)
}
