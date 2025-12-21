package zfactor

import (
	"math"
	"math/cmplx"
	"testing"
)

func almostEqualComplex(a, b complex128, tol float64) bool {
	return cmplx.Abs(a-b) < tol
}

func TestSolveCubic(t *testing.T) {
	tests := []struct {
		name       string
		a, b, c, d float64
		wantRoots  []complex128
		wantErr    bool
	}{
		{
			name: "Not cubic (a=0)",
			a:    0, b: 1, c: 2, d: 3,
			wantRoots: nil,
			wantErr:   true,
		},
		{
			name: "x^3 - 1 = 0",
			a:    1, b: 0, c: 0, d: -1,
			wantRoots: []complex128{1, complex(-0.5, math.Sqrt(3)/2), complex(-0.5, -math.Sqrt(3)/2)},
			wantErr:   false,
		},
		{
			name: "x^3 - 6x^2 + 11x - 6 = 0 (roots 1,2,3)",
			a:    1, b: -6, c: 11, d: -6,
			wantRoots: []complex128{1, 2, 3},
			wantErr:   false,
		},
		{
			name: "x^3 + 3x^2 + 3x + 1 = 0 (triple root -1)",
			a:    1, b: 3, c: 3, d: 1,
			wantRoots: []complex128{-1, -1, -1},
			wantErr:   false,
		},
	}

	const tol = 1e-6

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := SolveCubic(tt.a, tt.b, tt.c, tt.d)
			if tt.wantErr {
				if err == nil {
					t.Errorf("expected error but got nil")
				}
				return
			}
			if err != nil {
				t.Errorf("unexpected error: %v", err)
				return
			}

			for _, want := range tt.wantRoots {
				found := false
				for _, g := range got {
					if almostEqualComplex(g, want, tol) {
						found = true
						break
					}
				}
				if !found {
					t.Errorf("expected root %v not found in results %v", want, got)
				}
			}
		})
	}
}
