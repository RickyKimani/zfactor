package linalglite_test

import (
	"errors"
	"math"
	"math/rand"
	"testing"

	"github.com/rickykimani/zfactor/linalglite"
)

// randomMatrix returns a matrix of the given order with entries spanning
// several orders of magnitude, which is where elimination loses the most
// precision and so where a factorisation is worth testing.
func randomMatrix(rng *rand.Rand, n int) *linalglite.Dense {
	m := linalglite.New(n, n)

	for i := range n {
		row := m.Row(i)
		for j := range row {
			row[j] = (rng.Float64()*2 - 1) * math.Pow(10, rng.Float64()*6-3)
		}
	}

	return m
}

// randomVector returns a vector of the given length.
func randomVector(rng *rand.Rand, n int) []float64 {
	v := make([]float64, n)
	for i := range v {
		v[i] = (rng.Float64()*2 - 1) * math.Pow(10, rng.Float64()*4-2)
	}

	return v
}

// scale returns a figure to judge a residual against: the largest
// magnitude the multiplication had to work with.
func scale(a *linalglite.Dense, x, b []float64) float64 {
	largest := 1.0

	for _, v := range a.RawRowMajor() {
		if m := math.Abs(v); m > largest {
			largest = m
		}
	}

	var xLargest float64
	for _, v := range x {
		if m := math.Abs(v); m > xLargest {
			xLargest = m
		}
	}

	for _, v := range b {
		if m := math.Abs(v); m > largest {
			largest = m
		}
	}

	return largest * math.Max(xLargest, 1)
}

// TestSolveSatisfiesTheEquations is the central property: a returned
// solution, substituted back, must reproduce the right-hand side.
//
// It says nothing about how the answer was reached, which is the point.
// Every order is covered, so the closed forms for one, two and three are
// held to the same standard as the factorisation above them.
func TestSolveSatisfiesTheEquations(t *testing.T) {
	const relTol = 1e-10

	rng := rand.New(rand.NewSource(1))

	for _, n := range []int{1, 2, 3, 4, 5, 8, 16, 32, 64} {
		for range 200 {
			a := randomMatrix(rng, n)
			b := randomVector(rng, n)

			x, err := linalglite.Solve(a, b)
			if err != nil {
				if errors.Is(err, linalglite.ErrSingular) {
					continue
				}

				t.Fatalf("order %d: Solve returned an unexpected error: %v", n, err)
			}

			residual, err := linalglite.Residual(a, x, b)
			if err != nil {
				t.Fatalf("order %d: Residual returned an unexpected error: %v", n, err)
			}

			if rel := residual / scale(a, x, b); rel > relTol {
				t.Fatalf(
					"order %d: the solution leaves a residual of %.3e, %.3e relative",
					n, residual, rel,
				)
			}
		}
	}
}

// TestSolveAgainstAConstructedSolution checks against an answer known in
// advance rather than against the equations.
//
// A random matrix is multiplied by a chosen x to produce b, so the
// solution is known exactly before the solver sees it. Reproducing the
// equations and reproducing the intended answer are different claims, and
// an ill-conditioned matrix can satisfy the first while missing the
// second.
func TestSolveAgainstAConstructedSolution(t *testing.T) {
	const relTol = 1e-8

	rng := rand.New(rand.NewSource(2))

	for _, n := range []int{1, 2, 3, 5, 10, 25} {
		for range 100 {
			a := randomMatrix(rng, n)
			want := randomVector(rng, n)

			b, err := a.MulVec(want)
			if err != nil {
				t.Fatalf("MulVec returned an unexpected error: %v", err)
			}

			got, err := linalglite.Solve(a, b)
			if err != nil {
				if errors.Is(err, linalglite.ErrSingular) {
					continue
				}

				t.Fatalf("order %d: Solve returned an unexpected error: %v", n, err)
			}

			// Compare against the size of the intended answer, since a
			// well-conditioned system should recover it closely.
			var largest float64
			for _, v := range want {
				largest = math.Max(largest, math.Abs(v))
			}

			for i := range want {
				if rel := math.Abs(got[i]-want[i]) / math.Max(largest, 1); rel > relTol {
					// Ill-conditioning, not a defect, if the residual is
					// small: the system does not determine x closely.
					residual, _ := linalglite.Residual(a, got, b)
					if residual/scale(a, got, b) < 1e-12 {
						continue
					}

					t.Fatalf(
						"order %d, element %d: got %.9g, want %.9g",
						n, i, got[i], want[i],
					)
				}
			}
		}
	}
}

// TestInverseMultipliesBackToTheIdentity checks the inverse against its
// definition.
//
// Forming an inverse and multiplying by it is the least accurate way to
// use a factorisation, which is why Inverse is documented as a last
// resort and Solve preferred. The tolerance therefore scales with the
// conditioning of the matrix, estimated from the size of the inverse
// itself: a matrix that nearly fails to be invertible has a large
// inverse, and the product loses digits in proportion.
func TestInverseMultipliesBackToTheIdentity(t *testing.T) {
	rng := rand.New(rand.NewSource(3))

	for _, n := range []int{2, 3, 4, 7, 13, 32} {
		for range 20 {
			a := randomMatrix(rng, n)

			lu, err := linalglite.Factorize(a)
			if err != nil {
				if errors.Is(err, linalglite.ErrSingular) {
					continue
				}

				t.Fatalf("order %d: Factorize returned an unexpected error: %v", n, err)
			}

			inverse, err := lu.Inverse()
			if err != nil {
				t.Fatalf("order %d: Inverse returned an unexpected error: %v", n, err)
			}

			// An estimate of how much the product can lose: the largest
			// entry of the matrix times the largest of its inverse.
			var aLargest, iLargest float64
			for _, v := range a.RawRowMajor() {
				aLargest = math.Max(aLargest, math.Abs(v))
			}
			for _, v := range inverse.RawRowMajor() {
				iLargest = math.Max(iLargest, math.Abs(v))
			}

			tol := 1e-13 * aLargest * iLargest * float64(n)

			for j := range n {
				column := make([]float64, n)
				for i := range n {
					column[i] = inverse.At(i, j)
				}

				product, err := a.MulVec(column)
				if err != nil {
					t.Fatalf("MulVec returned an unexpected error: %v", err)
				}

				for i := range n {
					want := 0.0
					if i == j {
						want = 1
					}

					if math.Abs(product[i]-want) > tol {
						t.Fatalf(
							"order %d: A·A⁻¹ has %.3e at (%d, %d); want %g within %.3e",
							n, product[i], i, j, want, tol,
						)
					}
				}
			}
		}
	}
}

// TestDetAgainstKnownValues checks the determinant against matrices whose
// value is known by hand.
func TestDetAgainstKnownValues(t *testing.T) {
	const tol = 1e-9

	testCases := []struct {
		name string
		n    int
		data []float64
		want float64
	}{
		{"order one", 1, []float64{7}, 7},
		{"order two", 2, []float64{1, 2, 3, 4}, -2},
		{"order three", 3, []float64{6, 1, 1, 4, -2, 5, 2, 8, 7}, -306},
		{"identity of order four", 4, []float64{
			1, 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1,
		}, 1},
		{"upper triangular of order four", 4, []float64{
			2, 9, 9, 9,
			0, 3, 9, 9,
			0, 0, 4, 9,
			0, 0, 0, 5,
		}, 120},
		{"a swapped pair changes the sign", 2, []float64{3, 4, 1, 2}, 2},
		{"a repeated row is singular", 3, []float64{
			1, 2, 3,
			1, 2, 3,
			4, 5, 6,
		}, 0},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			a := linalglite.NewFrom(tc.n, tc.n, tc.data)

			got, err := linalglite.Det(a)
			if err != nil {
				t.Fatalf("Det returned an unexpected error: %v", err)
			}

			if math.Abs(got-tc.want) > tol*math.Max(math.Abs(tc.want), 1) {
				t.Errorf("determinant = %.9g; want %g", got, tc.want)
			}
		})
	}
}

// TestDetMatchesTheProductOfEigenvaluesOfATriangle checks the determinant
// of a triangular matrix, whose value is the product of its diagonal
// however large the matrix.
//
// It is a case where the answer is known for any order, so it tests the
// factorisation path rather than the closed forms.
func TestDetMatchesTheDiagonalOfATriangle(t *testing.T) {
	const relTol = 1e-12

	rng := rand.New(rand.NewSource(4))

	for _, n := range []int{4, 6, 10, 20} {
		a := linalglite.New(n, n)

		want := 1.0

		for i := range n {
			row := a.Row(i)

			// Diagonal entries away from zero, so the product is well
			// determined.
			d := rng.Float64() + 0.5
			row[i] = d
			want *= d

			for j := i + 1; j < n; j++ {
				row[j] = rng.Float64()*2 - 1
			}
		}

		got, err := linalglite.Det(a)
		if err != nil {
			t.Fatalf("order %d: Det returned an unexpected error: %v", n, err)
		}

		if rel := math.Abs(got-want) / math.Abs(want); rel > relTol {
			t.Errorf("order %d: determinant = %.12g; want %.12g", n, got, want)
		}
	}
}

// TestLogAbsDetSurvivesOverflow checks the logarithmic form on a matrix
// whose determinant is not representable.
//
// The product of the diagonal overflows to infinity while its logarithm
// remains an ordinary number, which is the reason the second form exists.
func TestLogAbsDetSurvivesOverflow(t *testing.T) {
	const n = 40

	a := linalglite.New(n, n)
	for i := range n {
		// Each diagonal entry large enough that the product exceeds the
		// range of a float64.
		a.Set(i, i, 1e20)
	}

	lu, err := linalglite.Factorize(a)
	if err != nil {
		t.Fatalf("Factorize returned an unexpected error: %v", err)
	}

	if det := lu.Det(); !math.IsInf(det, 1) {
		t.Errorf("the direct product gives %g; expected it to overflow", det)
	}

	log, sign, want := 0.0, 0.0, float64(n)*math.Log(1e20)
	log, sign = lu.LogAbsDet()

	if sign != 1 {
		t.Errorf("sign = %g; want 1", sign)
	}

	if rel := math.Abs(log-want) / want; rel > 1e-12 {
		t.Errorf("log of the determinant = %.12g; want %.12g", log, want)
	}
}

// TestPartialPivotingKeepsAccuracy checks the reason the pivot search is
// there.
//
// The matrix below has a tiny leading entry. Eliminating with it divides
// the rest of the matrix by that small number, and without an interchange
// the precision lost is not recovered. Bringing the larger entry to the
// diagonal first keeps the solution accurate.
func TestPartialPivotingKeepsAccuracy(t *testing.T) {
	const relTol = 1e-12

	// x + y = 2 with a nearly singular first equation.
	a := linalglite.NewFrom(2, 2, []float64{
		1e-18, 1,
		1, 1,
	})
	b := []float64{1, 2}

	x, err := linalglite.Solve(a, b)
	if err != nil {
		t.Fatalf("Solve returned an unexpected error: %v", err)
	}

	residual, err := linalglite.Residual(a, x, b)
	if err != nil {
		t.Fatalf("Residual returned an unexpected error: %v", err)
	}

	if residual > relTol {
		t.Errorf("the solution leaves a residual of %.3e; pivoting should keep it near rounding", residual)
	}

	// Larger, and through the factorisation rather than the closed form.
	large := linalglite.NewFrom(4, 4, []float64{
		1e-20, 1, 2, 3,
		1, 2, 3, 4,
		2, 3, 5, 7,
		3, 5, 7, 11,
	})
	rhs := []float64{1, 2, 3, 4}

	x, err = linalglite.Solve(large, rhs)
	if err != nil {
		t.Fatalf("Solve returned an unexpected error: %v", err)
	}

	residual, err = linalglite.Residual(large, x, rhs)
	if err != nil {
		t.Fatalf("Residual returned an unexpected error: %v", err)
	}

	if rel := residual / scale(large, x, rhs); rel > 1e-12 {
		t.Errorf("the solution leaves a relative residual of %.3e", rel)
	}
}

// TestSingularMatrices checks that a system without a unique solution is
// reported rather than answered.
//
// Returning infinities from a vanishing pivot would produce values that
// look like an answer, which is worse than a refusal.
func TestSingularMatrices(t *testing.T) {
	testCases := []struct {
		name string
		n    int
		data []float64
	}{
		{"a zero element", 1, []float64{0}},
		{"a repeated row", 2, []float64{1, 2, 1, 2}},
		{"a proportional row", 2, []float64{1, 2, 3, 6}},
		{"a zero row", 3, []float64{1, 2, 3, 0, 0, 0, 4, 5, 6}},
		{"a zero column", 3, []float64{1, 0, 3, 4, 0, 6, 7, 0, 9}},
		{"a dependent third row", 3, []float64{1, 2, 3, 4, 5, 6, 5, 7, 9}},
		{"order four, rank three", 4, []float64{
			1, 2, 3, 4,
			2, 4, 6, 8,
			1, 0, 1, 0,
			0, 1, 0, 1,
		}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			a := linalglite.NewFrom(tc.n, tc.n, tc.data)
			b := make([]float64, tc.n)
			for i := range b {
				b[i] = 1
			}

			_, err := linalglite.Solve(a, b)

			if !errors.Is(err, linalglite.ErrSingular) {
				t.Errorf("Solve returned %v; want ErrSingular", err)
			}

			// The determinant of a singular matrix is a value, not an
			// error, and it is zero.
			det, err := linalglite.Det(a)
			if err != nil {
				t.Errorf("Det returned an unexpected error: %v", err)
			}

			if det != 0 {
				t.Errorf("determinant = %g; want 0", det)
			}
		})
	}
}

// TestSolveIntoAllocatesNothingAndPreservesInput checks the properties an
// iteration depends on.
//
// A Newton loop calls this once per step, so it must not allocate; and it
// must leave the right-hand side alone, since the caller often still needs
// it. Aliasing the destination onto the input is allowed and is the
// convenient form for overwriting a residual with a step.
func TestSolveIntoAllocatesNothingAndPreservesInput(t *testing.T) {
	a := linalglite.NewFrom(4, 4, []float64{
		4, 1, 0, 0,
		1, 4, 1, 0,
		0, 1, 4, 1,
		0, 0, 1, 4,
	})
	b := []float64{1, 2, 3, 4}

	lu, err := linalglite.Factorize(a)
	if err != nil {
		t.Fatalf("Factorize returned an unexpected error: %v", err)
	}

	original := make([]float64, len(b))
	copy(original, b)

	x := make([]float64, len(b))

	allocations := testing.AllocsPerRun(100, func() {
		if err := lu.SolveInto(x, b); err != nil {
			t.Fatalf("SolveInto returned an unexpected error: %v", err)
		}
	})

	if allocations != 0 {
		t.Errorf("SolveInto allocated %.1f times per call; want none", allocations)
	}

	for i := range b {
		if b[i] != original[i] {
			t.Errorf("the right-hand side was modified: %v, was %v", b, original)
			break
		}
	}

	// Solving in place must give the same answer.
	inPlace := make([]float64, len(b))
	copy(inPlace, b)

	if err := lu.SolveInto(inPlace, inPlace); err != nil {
		t.Fatalf("SolveInto in place returned an unexpected error: %v", err)
	}

	for i := range x {
		if math.Abs(inPlace[i]-x[i]) > 1e-15 {
			t.Errorf("solving in place gave %v; want %v", inPlace, x)
			break
		}
	}
}

// TestRefactorizeReusesStorage checks the path an iteration on a changing
// matrix takes.
//
// The matrix differs each step in a Newton method, so the factorisation
// must be replaceable without allocating afresh.
func TestRefactorizeReusesStorage(t *testing.T) {
	first := linalglite.NewFrom(3, 3, []float64{
		2, 0, 0,
		0, 3, 0,
		0, 0, 4,
	})

	lu, err := linalglite.Factorize(first)
	if err != nil {
		t.Fatalf("Factorize returned an unexpected error: %v", err)
	}

	if det := lu.Det(); math.Abs(det-24) > 1e-12 {
		t.Errorf("determinant = %g; want 24", det)
	}

	second := linalglite.NewFrom(3, 3, []float64{
		1, 2, 3,
		4, 5, 6,
		7, 8, 10,
	})

	allocations := testing.AllocsPerRun(100, func() {
		if err := lu.Refactorize(second); err != nil {
			t.Fatalf("Refactorize returned an unexpected error: %v", err)
		}
	})

	if allocations != 0 {
		t.Errorf("Refactorize allocated %.1f times per call; want none", allocations)
	}

	// The determinant of the second matrix is -3.
	if det := lu.Det(); math.Abs(det+3) > 1e-12 {
		t.Errorf("determinant after refactorising = %g; want -3", det)
	}

	// And it must solve the new system, not the old one.
	b := []float64{6, 15, 25}

	x, err := lu.Solve(b)
	if err != nil {
		t.Fatalf("Solve returned an unexpected error: %v", err)
	}

	residual, err := linalglite.Residual(second, x, b)
	if err != nil {
		t.Fatalf("Residual returned an unexpected error: %v", err)
	}

	if residual > 1e-10 {
		t.Errorf("the solution leaves a residual of %.3e", residual)
	}
}

// TestClosedFormsAgreeWithTheFactorisation checks the two paths against
// each other.
//
// Orders one to three are solved by Cramer's rule and larger ones by
// elimination. They are separate code, so agreement is the check that the
// shortcut is a shortcut and not a different answer.
func TestClosedFormsAgreeWithTheFactorisation(t *testing.T) {
	const relTol = 1e-9

	rng := rand.New(rand.NewSource(5))

	for _, n := range []int{1, 2, 3} {
		for range 500 {
			a := randomMatrix(rng, n)
			b := randomVector(rng, n)

			viaClosedForm, err := linalglite.Solve(a, b)
			if err != nil {
				continue
			}

			lu, err := linalglite.Factorize(a)
			if err != nil {
				continue
			}

			viaFactorisation, err := lu.Solve(b)
			if err != nil {
				t.Fatalf("order %d: Solve returned an unexpected error: %v", n, err)
			}

			var largest float64
			for _, v := range viaFactorisation {
				largest = math.Max(largest, math.Abs(v))
			}

			for i := range viaClosedForm {
				diff := math.Abs(viaClosedForm[i] - viaFactorisation[i])

				if rel := diff / math.Max(largest, 1); rel > relTol {
					t.Fatalf(
						"order %d, element %d: the closed form gives %.12g and the factorisation %.12g",
						n, i, viaClosedForm[i], viaFactorisation[i],
					)
				}
			}
		}
	}
}

// TestShapeErrors checks that mismatched sizes are refused.
func TestShapeErrors(t *testing.T) {
	square := linalglite.New(3, 3)
	oblong := linalglite.New(2, 3)

	t.Run("solving a non-square system", func(t *testing.T) {
		if _, err := linalglite.Solve(oblong, []float64{1, 2}); !errors.Is(err, linalglite.ErrShape) {
			t.Errorf("got %v; want ErrShape", err)
		}
	})

	t.Run("factorising a non-square matrix", func(t *testing.T) {
		if _, err := linalglite.Factorize(oblong); !errors.Is(err, linalglite.ErrShape) {
			t.Errorf("got %v; want ErrShape", err)
		}
	})

	t.Run("a determinant of a non-square matrix", func(t *testing.T) {
		if _, err := linalglite.Det(oblong); !errors.Is(err, linalglite.ErrShape) {
			t.Errorf("got %v; want ErrShape", err)
		}
	})

	t.Run("the wrong number of right-hand side values", func(t *testing.T) {
		if _, err := linalglite.Solve(square, []float64{1, 2}); !errors.Is(err, linalglite.ErrShape) {
			t.Errorf("got %v; want ErrShape", err)
		}
	})

	t.Run("multiplying by a vector of the wrong length", func(t *testing.T) {
		if _, err := square.MulVec([]float64{1, 2}); !errors.Is(err, linalglite.ErrShape) {
			t.Errorf("got %v; want ErrShape", err)
		}
	})

	t.Run("refactorising at a different order", func(t *testing.T) {
		lu, err := linalglite.Factorize(linalglite.Identity(3))
		if err != nil {
			t.Fatalf("Factorize returned an unexpected error: %v", err)
		}

		if err := lu.Refactorize(linalglite.Identity(4)); !errors.Is(err, linalglite.ErrShape) {
			t.Errorf("got %v; want ErrShape", err)
		}
	})

	t.Run("a destination of the wrong length", func(t *testing.T) {
		lu, err := linalglite.Factorize(square.Clone())
		if err != nil {
			// The zero matrix is singular; use an identity instead.
			lu, err = linalglite.Factorize(linalglite.Identity(3))
			if err != nil {
				t.Fatalf("Factorize returned an unexpected error: %v", err)
			}
		}

		if err := lu.SolveInto(make([]float64, 2), []float64{1, 2, 3}); !errors.Is(err, linalglite.ErrShape) {
			t.Errorf("got %v; want ErrShape", err)
		}
	})
}

// TestDenseAccessors checks the matrix type itself.
func TestDenseAccessors(t *testing.T) {
	m := linalglite.NewFrom(2, 3, []float64{1, 2, 3, 4, 5, 6})

	if m.Rows() != 2 || m.Cols() != 3 {
		t.Errorf("dimensions = %dx%d; want 2x3", m.Rows(), m.Cols())
	}

	if m.IsSquare() {
		t.Error("a 2x3 matrix reports itself square")
	}

	if got := m.At(1, 2); got != 6 {
		t.Errorf("At(1, 2) = %g; want 6", got)
	}

	m.Set(0, 0, 99)
	if got := m.At(0, 0); got != 99 {
		t.Errorf("At(0, 0) = %g after Set; want 99", got)
	}

	// A row is a view, so writing through it changes the matrix.
	row := m.Row(1)
	row[0] = 42

	if got := m.At(1, 0); got != 42 {
		t.Errorf("At(1, 0) = %g after writing through Row; want 42", got)
	}

	// NewFrom copies, so the caller's slice is independent.
	data := []float64{1, 2, 3, 4}
	n := linalglite.NewFrom(2, 2, data)
	data[0] = 99

	if got := n.At(0, 0); got != 1 {
		t.Errorf("At(0, 0) = %g; NewFrom should have copied", got)
	}

	// Clone is independent too.
	c := n.Clone()
	c.Set(0, 0, 7)

	if got := n.At(0, 0); got != 1 {
		t.Errorf("the original changed to %g when its clone was written", got)
	}

	c.Zero()
	for _, v := range c.RawRowMajor() {
		if v != 0 {
			t.Errorf("Zero left %g behind", v)
			break
		}
	}

	if got := linalglite.Identity(3); got.At(0, 0) != 1 || got.At(0, 1) != 0 || got.At(2, 2) != 1 {
		t.Error("Identity is not the identity")
	}
}

// TestDensePanics checks the cases documented as programming errors.
func TestDensePanics(t *testing.T) {
	testCases := []struct {
		name string
		call func()
	}{
		{"a non-positive dimension", func() { linalglite.New(0, 3) }},
		{"a negative dimension", func() { linalglite.New(-1, 3) }},
		{"the wrong number of elements", func() { linalglite.NewFrom(2, 2, []float64{1, 2, 3}) }},
		{"reading out of range", func() { linalglite.New(2, 2).At(2, 0) }},
		{"writing out of range", func() { linalglite.New(2, 2).Set(0, 2, 1) }},
		{"a row out of range", func() { linalglite.New(2, 2).Row(5) }},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			defer func() {
				if recover() == nil {
					t.Error("expected a panic; the call returned")
				}
			}()

			tc.call()
		})
	}
}

// TestMulVecIntoAllocatesNothing checks the path used to form a residual
// inside an iteration.
func TestMulVecIntoAllocatesNothing(t *testing.T) {
	a := linalglite.NewFrom(3, 3, []float64{1, 2, 3, 4, 5, 6, 7, 8, 10})
	x := []float64{1, 1, 1}
	dst := make([]float64, 3)

	allocations := testing.AllocsPerRun(100, func() {
		if err := a.MulVecInto(dst, x); err != nil {
			t.Fatalf("MulVecInto returned an unexpected error: %v", err)
		}
	})

	if allocations != 0 {
		t.Errorf("MulVecInto allocated %.1f times per call; want none", allocations)
	}

	for i, want := range []float64{6, 15, 25} {
		if dst[i] != want {
			t.Errorf("element %d = %g; want %g", i, dst[i], want)
		}
	}
}
