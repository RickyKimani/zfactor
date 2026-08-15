package linalglite

import (
	"errors"
	"math"
	"math/rand"
	"testing"
)

// TestDecompositionReconstructsTheMatrix checks the factorisation against
// its defining property, reading the factors directly.
//
// Partial pivoting produces P·A = L·U, so applying the recorded
// interchanges to the original matrix must reproduce the product of the
// two triangles. This is the direct check: it reaches the stored factors
// and the pivot record rather than inferring their correctness from
// something solved with them, and it would catch an error in either that
// a solve could mask.
//
// It lives inside the package because the factors are not exported, and
// they are not exported because nothing outside needs them.
func TestDecompositionReconstructsTheMatrix(t *testing.T) {
	const relTol = 1e-13

	rng := rand.New(rand.NewSource(11))

	for _, n := range []int{2, 3, 4, 7, 13, 32, 50} {
		for range 50 {
			a := New(n, n)
			for i := range n {
				row := a.Row(i)
				for j := range row {
					row[j] = (rng.Float64()*2 - 1) * math.Pow(10, rng.Float64()*4-2)
				}
			}

			original := a.Clone()

			lu, err := Factorize(a)
			if err != nil {
				if errors.Is(err, ErrSingular) {
					continue
				}

				t.Fatalf("order %d: Factorize returned an unexpected error: %v", n, err)
			}

			// Apply the interchanges to a copy of the original, in the
			// order they were made, giving P·A.
			permuted := original.Clone()
			for k := range n {
				if p := lu.pivot[k]; p != k {
					rowK := permuted.Row(k)
					rowP := permuted.Row(p)

					for j := range rowK {
						rowK[j], rowP[j] = rowP[j], rowK[j]
					}
				}
			}

			// Form L·U from the stored factors. L has an implicit unit
			// diagonal below which the multipliers sit; U is on and above
			// the diagonal.
			var largest float64
			for _, v := range original.RawRowMajor() {
				largest = math.Max(largest, math.Abs(v))
			}

			for i := range n {
				for j := range n {
					var sum float64

					// L[i][k] is one at k == i, the stored value below it,
					// and zero above; U[k][j] is zero for k > j.
					for k := 0; k <= min(i, j); k++ {
						l := 1.0
						if k < i {
							l = lu.lu[i*n+k]
						}

						sum += l * lu.lu[k*n+j]
					}

					want := permuted.At(i, j)

					if diff := math.Abs(sum - want); diff > relTol*math.Max(largest, 1)*float64(n) {
						t.Fatalf(
							"order %d at (%d, %d): L·U gives %.12g but P·A is %.12g",
							n, i, j, sum, want,
						)
					}
				}
			}
		}
	}
}

// TestDecompositionLeavesTheInputAlone checks that factorising copies
// rather than consuming the matrix.
//
// A caller forming a Jacobian, factorising it and then wanting to compute
// a residual against it needs the matrix intact.
func TestDecompositionLeavesTheInputAlone(t *testing.T) {
	a := NewFrom(3, 3, []float64{
		1e-15, 2, 3,
		4, 5, 6,
		7, 8, 10,
	})

	before := a.Clone()

	if _, err := Factorize(a); err != nil {
		t.Fatalf("Factorize returned an unexpected error: %v", err)
	}

	for i, v := range a.RawRowMajor() {
		if v != before.RawRowMajor()[i] {
			t.Fatalf("the matrix was modified: %v, was %v", a.RawRowMajor(), before.RawRowMajor())
		}
	}
}

// TestPivotsSelectTheLargestEntry checks that the interchange chooses the
// largest available pivot rather than merely a non-zero one.
//
// The magnitude of the pivot is what governs how much precision the
// elimination loses, so picking the first usable entry would be
// materially worse than picking the biggest.
func TestPivotsSelectTheLargestEntry(t *testing.T) {
	// The largest entry of the first column is in the last row.
	a := NewFrom(3, 3, []float64{
		1, 1, 1,
		2, 5, 3,
		9, 2, 4,
	})

	lu, err := Factorize(a)
	if err != nil {
		t.Fatalf("Factorize returned an unexpected error: %v", err)
	}

	if lu.pivot[0] != 2 {
		t.Errorf("the first pivot is row %d; the largest entry is in row 2", lu.pivot[0])
	}

	// Every multiplier stored below the diagonal must be no larger than
	// one, which is what pivoting on the largest entry guarantees and the
	// reason the factorisation stays stable.
	for i := range lu.n {
		for k := range i {
			if m := math.Abs(lu.lu[i*lu.n+k]); m > 1+1e-15 {
				t.Errorf("the multiplier at (%d, %d) is %g; pivoting should keep it within one", i, k, m)
			}
		}
	}
}

// TestMultipliersStayBoundedOnRandomMatrices extends that guarantee
// across many matrices.
//
// It is the property that distinguishes partial pivoting from none: every
// multiplier is a ratio of an entry to the largest in its column, so none
// exceeds one, and the elimination cannot amplify what it subtracts.
func TestMultipliersStayBoundedOnRandomMatrices(t *testing.T) {
	rng := rand.New(rand.NewSource(12))

	for _, n := range []int{4, 10, 25} {
		for range 100 {
			a := New(n, n)
			for i := range n {
				row := a.Row(i)
				for j := range row {
					row[j] = (rng.Float64()*2 - 1) * math.Pow(10, rng.Float64()*6-3)
				}
			}

			lu, err := Factorize(a)
			if err != nil {
				continue
			}

			for i := range n {
				for k := range i {
					if m := math.Abs(lu.lu[i*n+k]); m > 1+1e-12 {
						t.Fatalf("order %d: the multiplier at (%d, %d) is %g", n, i, k, m)
					}
				}
			}
		}
	}
}

// TestSparseUpdateSkipIsSoundlyTaken checks the shortcut in the inner
// loop.
//
// A zero multiplier means the row needs no update, and the loop skips it.
// Stoichiometric and Jacobian matrices in this field are often sparse
// enough for that to matter, so the result must be the same as if the
// subtraction had been carried out.
func TestSparseUpdateSkipIsSoundlyTaken(t *testing.T) {
	// A matrix with structural zeros below the diagonal, so the skip is
	// taken repeatedly.
	a := NewFrom(5, 5, []float64{
		4, 1, 0, 0, 0,
		0, 4, 1, 0, 0,
		0, 0, 4, 1, 0,
		0, 0, 0, 4, 1,
		1, 0, 0, 0, 4,
	})

	b := []float64{1, 2, 3, 4, 5}

	x, err := Solve(a, b)
	if err != nil {
		t.Fatalf("Solve returned an unexpected error: %v", err)
	}

	residual, err := Residual(a, x, b)
	if err != nil {
		t.Fatalf("Residual returned an unexpected error: %v", err)
	}

	if residual > 1e-13 {
		t.Errorf("the solution leaves a residual of %.3e", residual)
	}
}
