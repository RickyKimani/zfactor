package zfactor_test

import (
	"errors"
	"fmt"
	"math/cmplx"
	"strings"
	"testing"

	"github.com/rickykimani/zfactor"
)

// TestRangeErrorMessage checks both forms of the message.
//
// The type is shared by correlations that know which substance they
// describe and by those that do not, so the name is omitted rather than
// rendered empty when it is absent.
func TestRangeErrorMessage(t *testing.T) {
	t.Run("with a name", func(t *testing.T) {
		err := &zfactor.RangeError{
			Name: "Carbon dioxide",
			T:    293.15,
			Low:  298.15,
			High: 2000,
		}

		got := err.Error()

		for _, want := range []string{"293.15", "298.15", "2000", "Carbon dioxide"} {
			if !strings.Contains(got, want) {
				t.Errorf("message %q does not mention %q", got, want)
			}
		}
	})

	t.Run("without a name", func(t *testing.T) {
		err := &zfactor.RangeError{T: 250, Low: 300, High: 400}

		got := err.Error()

		if strings.Contains(got, "fitted for") {
			t.Errorf("message %q names a substance it was not given", got)
		}

		for _, want := range []string{"250", "300", "400"} {
			if !strings.Contains(got, want) {
				t.Errorf("message %q does not mention %q", got, want)
			}
		}
	})
}

// TestRangeErrorIsDiscoverable checks that the shared type can be
// recovered from a wrapped error, which is how callers distinguish an
// extrapolation from a failure.
func TestRangeErrorIsDiscoverable(t *testing.T) {
	wrapped := fmt.Errorf("computing a property: %w", &zfactor.RangeError{
		Name: "Water", T: 250, Low: 273.15, High: 373.15,
	})

	var rangeErr *zfactor.RangeError
	if !errors.As(wrapped, &rangeErr) {
		t.Fatalf("errors.As did not recover the range error from %v", wrapped)
	}

	if rangeErr.Name != "Water" || rangeErr.T != 250 {
		t.Errorf("recovered %+v; want the original values", rangeErr)
	}
}

// TestSolveCubicRepeatedRoots checks the polynomials whose derivative
// vanishes at a root.
//
// Newton's method cannot refine a repeated root, since the derivative it
// divides by is zero there, and the refinement is abandoned rather than
// producing an infinity. The roots Cardano supplies must therefore
// survive unchanged.
func TestSolveCubicRepeatedRoots(t *testing.T) {
	testCases := []struct {
		name       string
		a, b, c, d float64
		want       []complex128
	}{
		{"triple root at -1", 1, 3, 3, 1, []complex128{-1, -1, -1}},
		{"triple root at 0", 1, 0, 0, 0, []complex128{0, 0, 0}},
		{"double root at 1 with a simple root at -2", 1, 0, -3, 2, []complex128{1, 1, -2}},
	}

	const tol = 1e-6

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := zfactor.SolveCubic(tc.a, tc.b, tc.c, tc.d)
			if err != nil {
				t.Fatalf("SolveCubic returned an unexpected error: %v", err)
			}

			for _, want := range tc.want {
				var found bool
				for _, g := range got {
					if cmplx.Abs(g-want) < tol {
						found = true
						break
					}
				}

				if !found {
					t.Errorf("root %v is missing from %v", want, got)
				}
			}
		})
	}
}

// TestSolveCubicRefinesTheRoots checks that the refinement is applied,
// by comparing against the accuracy the closed form alone reaches.
//
// The polynomial below has roots spread over six orders of magnitude,
// which is where the substitution removing the quadratic term loses the
// most precision. Cardano's estimate leaves a residual far above
// rounding; refined, it should sit near it.
func TestSolveCubicRefinesTheRoots(t *testing.T) {
	const (
		a, b, c, d = 1.0, -1000001.0, 1000000.0, 0.0
		tol        = 1e-9
	)

	roots, err := zfactor.SolveCubic(a, b, c, d)
	if err != nil {
		t.Fatalf("SolveCubic returned an unexpected error: %v", err)
	}

	// The roots are 0, 1 and 1000000.
	for _, want := range []complex128{0, 1, 1000000} {
		var closest float64 = 1e300

		for _, got := range roots {
			scale := cmplx.Abs(want)
			if scale < 1 {
				scale = 1
			}

			if rel := cmplx.Abs(got-want) / scale; rel < closest {
				closest = rel
			}
		}

		if closest > tol {
			t.Errorf("root %v is not reproduced; nearest is %.3e away relatively (roots %v)", want, closest, roots)
		}
	}
}
