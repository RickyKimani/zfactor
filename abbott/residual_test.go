package abbott

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
)

// Carbon dioxide, from Table B.1 of Smith, Van Ness & Abbott.
const (
	co2Tc = 304.2
	co2W  = 0.224
)

// TestExample6_5FinalState reproduces step 3 of Example 6.5 of Smith,
// Van Ness & Abbott: the residual enthalpy and entropy of carbon dioxide
// at 20°C and 15 bar, obtained from the generalized second-virial-
// coefficient correlation.
//
// The example applies this correlation at the final state rather than
// the Lee-Kesler tables because the reduced pressure is low, where the
// truncated virial equation is adequate.
//
// The functions return dimensionless residuals, so both are scaled to
// J·mol⁻¹ and J·mol⁻¹·K⁻¹ for comparison with the published values.
func TestExample6_5FinalState(t *testing.T) {
	const (
		Tr = 0.964
		Pr = 0.203

		wantEnthalpy = -660.0 // J/mol
		wantEntropy  = -1.59  // J/(mol*K)

		// A relative tolerance, because the example quotes the reduced
		// conditions rounded to three decimals but carries more
		// precision through the calculation. Both residuals are linear
		// in Pr, so rounding it from 0.2032 to 0.203 alone accounts for
		// most of the 0.2% difference seen here.
		relTol = 5e-3
	)

	t.Run("residual enthalpy", func(t *testing.T) {
		got, err := ResidualEnthalpy(Tr, Pr, co2W)
		if err != nil {
			t.Fatalf("ResidualEnthalpy returned an unexpected error: %v", err)
		}

		h := zfactor.RSI * co2Tc * got

		if rel := math.Abs(h-wantEnthalpy) / math.Abs(wantEnthalpy); rel > relTol {
			t.Errorf("residual enthalpy = %.1f J/mol; want %.0f (%.2f%% apart)", h, wantEnthalpy, 100*rel)
		}
	})

	t.Run("residual entropy", func(t *testing.T) {
		got, err := ResidualEntropy(Tr, Pr, co2W)
		if err != nil {
			t.Fatalf("ResidualEntropy returned an unexpected error: %v", err)
		}

		s := zfactor.RSI * got

		if rel := math.Abs(s-wantEntropy) / math.Abs(wantEntropy); rel > relTol {
			t.Errorf("residual entropy = %.3f J/(mol*K); want %.2f (%.2f%% apart)", s, wantEntropy, 100*rel)
		}
	})
}

// TestDerivativeCorrelations checks the tabulated derivative correlations
// against numerical derivatives of the correlations they differentiate.
//
// DB0 and DB1 are published as separate closed forms rather than being
// derived in code, so nothing otherwise ties them to B0 and B1. An error
// in either would be invisible: the residual properties would simply be
// wrong.
//
// The agreement is limited by the published coefficients themselves,
// which are rounded to three significant figures. DB0 carries 0.675
// where differentiating B0 gives 0.422 × 1.6 = 0.6752, and DB1 carries
// 0.722 against 0.172 × 4.2 = 0.7224. Those account for relative offsets
// of 3.0e-4 and 5.5e-4, which are constant in Tr and set the tolerance
// here.
func TestDerivativeCorrelations(t *testing.T) {
	const (
		h   = 1e-7
		tol = 1e-3
	)

	testCases := []struct {
		name       string
		value      func(float64) (float64, error)
		derivative func(float64) (float64, error)
	}{
		{"B0", B0, DB0},
		{"B1", B1, DB1},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			for _, Tr := range []float64{0.5, 0.7, 0.964, 1.128, 1.5, 2.0, 3.0} {
				upper, err := tc.value(Tr + h)
				if err != nil {
					t.Fatalf("evaluation failed: %v", err)
				}

				lower, err := tc.value(Tr - h)
				if err != nil {
					t.Fatalf("evaluation failed: %v", err)
				}

				want := (upper - lower) / (2 * h)

				got, err := tc.derivative(Tr)
				if err != nil {
					t.Fatalf("derivative failed: %v", err)
				}

				if rel := math.Abs(got-want) / math.Abs(want); rel > tol {
					t.Errorf(
						"at Tr = %.3f: derivative = %.8f but the numerical derivative is %.8f (%.2e apart)",
						Tr, got, want, rel,
					)
				}
			}
		})
	}
}

// TestResidualsAreLinearInPressure checks that both residual properties
// scale in proportion to the reduced pressure.
//
// The truncated virial equation is linear in pressure, so the residual
// properties derived from it inherit that exactly. The relation holds
// regardless of temperature or acentric factor.
func TestResidualsAreLinearInPressure(t *testing.T) {
	const (
		Tr  = 0.964
		Pr  = 0.203
		tol = 1e-12
	)

	testCases := []struct {
		name string
		fn   func(Tr, Pr, acentric float64) (float64, error)
	}{
		{"residual enthalpy", ResidualEnthalpy},
		{"residual entropy", ResidualEntropy},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			single, err := tc.fn(Tr, Pr, co2W)
			if err != nil {
				t.Fatalf("evaluation failed: %v", err)
			}

			double, err := tc.fn(Tr, 2*Pr, co2W)
			if err != nil {
				t.Fatalf("evaluation failed: %v", err)
			}

			if math.Abs(double-2*single) > tol {
				t.Errorf(
					"doubling the reduced pressure gave %.12f; want twice %.12f",
					double, single,
				)
			}
		})
	}
}

// TestResidualsAtZeroAcentricFactor checks that a simple fluid reduces to
// the base term alone, since the departure term is weighted by the
// acentric factor.
func TestResidualsAtZeroAcentricFactor(t *testing.T) {
	const (
		Tr  = 0.964
		Pr  = 0.203
		tol = 1e-12
	)

	b0, err := B0(Tr)
	if err != nil {
		t.Fatalf("B0 failed: %v", err)
	}

	db0, err := DB0(Tr)
	if err != nil {
		t.Fatalf("DB0 failed: %v", err)
	}

	t.Run("residual enthalpy", func(t *testing.T) {
		got, err := ResidualEnthalpy(Tr, Pr, 0)
		if err != nil {
			t.Fatalf("ResidualEnthalpy returned an unexpected error: %v", err)
		}

		if want := Pr * (b0 - Tr*db0); math.Abs(got-want) > tol {
			t.Errorf("got %.12f; want Pr*(B0 - Tr*DB0) = %.12f", got, want)
		}
	})

	t.Run("residual entropy", func(t *testing.T) {
		got, err := ResidualEntropy(Tr, Pr, 0)
		if err != nil {
			t.Fatalf("ResidualEntropy returned an unexpected error: %v", err)
		}

		if want := -Pr * db0; math.Abs(got-want) > tol {
			t.Errorf("got %.12f; want -Pr*DB0 = %.12f", got, want)
		}
	})
}

// TestResidualInvalidInput checks that non-physical reduced conditions
// are rejected.
func TestResidualInvalidInput(t *testing.T) {
	testCases := []struct {
		name   string
		Tr, Pr float64
	}{
		{"non-positive reduced temperature", 0, 0.203},
		{"negative reduced temperature", -1, 0.203},
		{"non-positive reduced pressure", 0.964, 0},
		{"negative reduced pressure", 0.964, -1},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := ResidualEnthalpy(tc.Tr, tc.Pr, co2W); err == nil {
				t.Error("ResidualEnthalpy: expected an error; got nil")
			}

			if _, err := ResidualEntropy(tc.Tr, tc.Pr, co2W); err == nil {
				t.Error("ResidualEntropy: expected an error; got nil")
			}
		})
	}
}
