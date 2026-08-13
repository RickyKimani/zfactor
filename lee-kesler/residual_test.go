package leekesler

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

// TestExample6_5InitialState reproduces step 1 of Example 6.5 of Smith,
// Van Ness & Abbott: the residual enthalpy and entropy of supercritical
// carbon dioxide at 70°C and 150 bar.
//
// The example reads all four residual tables at Tr = 1.128 and
// Pr = 2.032 and combines each pair with the acentric factor. The
// Lee-Kesler tables are used rather than the virial correlation because
// the reduced pressure is high, where the truncated virial equation is
// no longer adequate.
//
// The enthalpy tables are dimensionless in H^R/(R·Tc) and the entropy
// tables in S^R/R, so each is scaled accordingly. The published step
// values are the negated residuals, since step 1 transforms the real
// fluid into its ideal-gas state.
func TestExample6_5InitialState(t *testing.T) {
	const (
		Tr = 1.128
		Pr = 2.032

		wantH0 = -2.709
		wantH1 = -0.921
		wantS0 = -1.846
		wantS1 = -0.938

		tol = 1e-3
	)

	t.Run("table terms", func(t *testing.T) {
		testCases := []struct {
			name  string
			table *table
			want  float64
		}{
			{"base residual enthalpy", H0Table, wantH0},
			{"departure residual enthalpy", H1Table, wantH1},
			{"base residual entropy", S0Table, wantS0},
			{"departure residual entropy", S1Table, wantS1},
		}

		for _, tc := range testCases {
			t.Run(tc.name, func(t *testing.T) {
				got, err := tc.table.At(Tr, Pr)
				if err != nil {
					t.Fatalf("At returned an unexpected error: %v", err)
				}

				if math.Abs(got-tc.want) > tol {
					t.Errorf("got %.6f; want %.3f", got, tc.want)
				}
			})
		}
	})

	// The step values are compared with a relative tolerance rather than
	// an absolute one. The example quotes the reduced conditions rounded
	// to three decimals but carries more precision through the
	// calculation, so reproducing it from the printed Tr and Pr shifts
	// the result slightly — by 0.03% here.
	const relTol = 1e-3

	t.Run("enthalpy change of step 1", func(t *testing.T) {
		const want = 7372.0 // J/mol

		residual, err := ResidualEnthalpy.Eval(Tr, Pr, co2W)
		if err != nil {
			t.Fatalf("Eval returned an unexpected error: %v", err)
		}

		got := -zfactor.RSI * co2Tc * residual

		if rel := math.Abs(got-want) / math.Abs(want); rel > relTol {
			t.Errorf("enthalpy change = %.1f J/mol; want %.0f (%.3f%% apart)", got, want, 100*rel)
		}
	})

	t.Run("entropy change of step 1", func(t *testing.T) {
		const want = 17.09 // J/(mol*K)

		residual, err := ResidualEntropy.Eval(Tr, Pr, co2W)
		if err != nil {
			t.Fatalf("Eval returned an unexpected error: %v", err)
		}

		got := -zfactor.RSI * residual

		if rel := math.Abs(got-want) / math.Abs(want); rel > relTol {
			t.Errorf("entropy change = %.3f J/(mol*K); want %.2f (%.3f%% apart)", got, want, 100*rel)
		}
	})
}

// TestResidualsVanishAtLowPressure checks that both residual properties
// approach zero as the pressure falls, since a gas approaches its
// ideal-gas state there and the residuals measure the departure from it.
//
// Only supercritical temperatures are checked, since below the critical
// temperature the low-pressure entries describe a condensed phase, whose
// residual properties remain large.
func TestResidualsVanishAtLowPressure(t *testing.T) {
	const tol = 0.05

	lowest := H0Table.Pr[0]

	for _, tr := range H0Table.Tr {
		if tr < 1 {
			continue
		}

		enthalpy, err := ResidualEnthalpy.Eval(tr, lowest, co2W)
		if err != nil {
			t.Fatalf("Eval returned an unexpected error: %v", err)
		}

		if math.Abs(enthalpy) > tol {
			t.Errorf(
				"at Tr = %g and Pr = %g: residual enthalpy = %.4f; want approximately 0",
				tr, lowest, enthalpy,
			)
		}

		entropy, err := ResidualEntropy.Eval(tr, lowest, co2W)
		if err != nil {
			t.Fatalf("Eval returned an unexpected error: %v", err)
		}

		if math.Abs(entropy) > tol {
			t.Errorf(
				"at Tr = %g and Pr = %g: residual entropy = %.4f; want approximately 0",
				tr, lowest, entropy,
			)
		}
	}
}

// TestResidualsAreNegativeForRealGases checks the sign of the residual
// properties across the supercritical region.
//
// Attractive intermolecular forces lower both the enthalpy and the
// entropy of a real gas relative to its ideal-gas state, so both
// residuals are negative wherever those forces dominate. A sign error in
// the combination would show here even though the magnitudes remain
// plausible.
func TestResidualsAreNegativeForRealGases(t *testing.T) {
	for _, tr := range []float64{1.0, 1.2, 1.5, 2.0} {
		for _, pr := range []float64{0.5, 1.0, 2.0, 3.0} {
			enthalpy, err := ResidualEnthalpy.Eval(tr, pr, co2W)
			if err != nil {
				t.Fatalf("Eval returned an unexpected error: %v", err)
			}

			if enthalpy >= 0 {
				t.Errorf(
					"at Tr = %.1f and Pr = %.1f: residual enthalpy = %.4f; want a negative value",
					tr, pr, enthalpy,
				)
			}

			entropy, err := ResidualEntropy.Eval(tr, pr, co2W)
			if err != nil {
				t.Fatalf("Eval returned an unexpected error: %v", err)
			}

			if entropy >= 0 {
				t.Errorf(
					"at Tr = %.1f and Pr = %.1f: residual entropy = %.4f; want a negative value",
					tr, pr, entropy,
				)
			}
		}
	}
}
