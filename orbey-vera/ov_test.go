package ov_test

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/abbott"
	ov "github.com/rickykimani/zfactor/orbey-vera"
	"github.com/rickykimani/zfactor/substance"
)

// TestC0 checks the simple-fluid contribution to the reduced third virial
// coefficient against Example 3.13 for n-butane at 470 K. It also verifies
// that non-positive reduced temperatures are rejected.
func TestC0(t *testing.T) {
	tr := 470 / substance.NButane.Critical.Tc

	tests := []struct {
		name    string
		tr      float64
		want    float64
		wantErr error
	}{
		{"Example 3.13 n-butane at 470 K", tr, 0.03498, nil},
		{"Invalid Tr=0", 0, 0, zfactor.ErrInvalidTr},
		{"Invalid Tr=-1", -1, 0, zfactor.ErrInvalidTr},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ov.C0(tt.tr)
			if err != tt.wantErr {
				t.Errorf("C0() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err == nil && math.Abs(got-tt.want)/tt.want > 5e-3 {
				t.Errorf("C0() = %.6f, want %.6f", got, tt.want)
			}
		})
	}
}

// TestC1 checks the acentric-factor contribution to the reduced third virial
// coefficient against Example 3.13 for n-butane at 470 K. It also verifies
// that non-positive reduced temperatures are rejected.
func TestC1(t *testing.T) {
	tr := 470 / substance.NButane.Critical.Tc

	tests := []struct {
		name    string
		tr      float64
		want    float64
		wantErr error
	}{
		{"Example 3.13 n-butane at 470 K", tr, 0.013724, nil},
		{"Invalid Tr=0", 0, 0, zfactor.ErrInvalidTr},
		{"Invalid Tr=-1", -1, 0, zfactor.ErrInvalidTr},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ov.C1(tt.tr)
			if err != tt.wantErr {
				t.Errorf("C1() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err == nil && math.Abs(got-tt.want)/tt.want > 5e-3 {
				t.Errorf("C1() = %.6f, want %.6f", got, tt.want)
			}
		})
	}
}

// TestReducedC checks the published reduced third virial coefficient from
// Example 3.13 and verifies its defining corresponding-states identity:
//
//	C^ = C0 + omega*C1
//
// The identity is checked independently of the rounded textbook value.
func TestReducedC(t *testing.T) {
	const T = 470.0
	nButane := substance.NButane
	tr := T / nButane.Critical.Tc
	w := nButane.Acentric

	c0, err := ov.C0(tr)
	if err != nil {
		t.Fatalf("C0 returned an unexpected error: %v", err)
	}
	c1, err := ov.C1(tr)
	if err != nil {
		t.Fatalf("C1 returned an unexpected error: %v", err)
	}

	got, err := ov.ReducedC(tr, w)
	if err != nil {
		t.Fatalf("ReducedC returned an unexpected error: %v", err)
	}

	if want := 0.03772; math.Abs(got-want)/want > 5e-3 {
		t.Errorf("ReducedC() = %.6f, want %.6f", got, want)
	}

	if want := c0 + w*c1; math.Abs(got-want) > 1e-12 {
		t.Errorf("ReducedC() = %.12f, want C0 + w*C1 = %.12f", got, want)
	}
}

// TestExample3_10Third reproduces the third-virial portion of Example 3.10
// for n-butane at 510 K and 25 bar. It checks the reduced second and third
// virial coefficients and the compressibility factor obtained from the
// smallest positive reduced-density root. The textbook values are rounded,
// so the comparisons use a relative tolerance.
func TestExample3_10Third(t *testing.T) {
	nButane := substance.NButane
	const (
		T = 510 //K
		P = 25  //bar

		wantB = -0.220
		wantC = 0.0352
		wantZ = 0.876

		relTol = 1e-2
	)

	Tr := T / nButane.Critical.Tc
	Pr := P / nButane.Critical.Pc
	w := nButane.Acentric

	c, err := ov.ReducedC(Tr, w)
	if err != nil {
		t.Fatalf("ReducedC returned an unexpected error: %v", err)
	}
	b, err := abbott.ReducedB(Tr, w)
	if err != nil {
		t.Fatalf("ReducedB returned an unexpected error: %v", err)
	}

	delta, err := ov.Delta(
		ov.DeltaArgs{
			C:  c,
			B:  b,
			Tr: Tr,
			Pr: Pr,
		},
	)
	if err != nil {
		t.Fatalf("Delta returned an unexpected error: %v", err)
	}
	z := Pr / (delta * Tr)

	checks := []struct {
		name string
		got  float64
		want float64
	}{
		{"reduced second virial coefficient", b, wantB},
		{"reduced third virial coefficient", c, wantC},
		{"compressibility factor", z, wantZ},
	}

	for _, check := range checks {
		t.Run(check.name, func(t *testing.T) {
			if rel := math.Abs(check.got-check.want) / math.Abs(check.want); rel > relTol {
				t.Errorf("got %.6f; want %.6f (%.2f%% apart)", check.got, check.want, 100*rel)
			}
		})
	}
}
