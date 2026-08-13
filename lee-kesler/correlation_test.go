package leekesler

import (
	"math"
	"testing"
)

// tables returns every generated correlation table with a label, so the
// structural invariants can be checked against all of them.
func tables() map[string]*table {
	return map[string]*table{
		"Z0":   Z0Table,
		"Z1":   Z1Table,
		"H0":   H0Table,
		"H1":   H1Table,
		"S0":   S0Table,
		"S1":   S1Table,
		"PHI0": PHI0Table,
		"PHI1": PHI1Table,
	}
}

// TestExample3_10 reproduces part (b) of Example 3.10 of Smith, Van Ness
// & Abbott: the molar volume of n-butane at 510 K and 25 bar from the
// generalized compressibility-factor correlation.
//
// The example works in reduced conditions Tr = 1.200 and Pr = 0.659,
// reads Z⁰ and Z¹ from the Lee-Kesler tables and combines them with the
// acentric factor of 0.200.
//
// The published intermediate values are quoted to three decimal places
// and rounded before being combined, so the tolerances here reflect that
// rounding rather than the accuracy of the interpolation.
func TestExample3_10(t *testing.T) {
	const (
		Tr = 1.200
		Pr = 0.659
		w  = 0.200

		wantZ0 = 0.865
		wantZ1 = 0.038
		wantZ  = 0.873

		tol = 1e-3
	)

	t.Run("base term", func(t *testing.T) {
		got, err := Z0Table.At(Tr, Pr)
		if err != nil {
			t.Fatalf("At returned an unexpected error: %v", err)
		}

		if math.Abs(got-wantZ0) > tol {
			t.Errorf("Z0 = %.6f; want %.3f", got, wantZ0)
		}
	})

	t.Run("departure term", func(t *testing.T) {
		got, err := Z1Table.At(Tr, Pr)
		if err != nil {
			t.Fatalf("At returned an unexpected error: %v", err)
		}

		if math.Abs(got-wantZ1) > tol {
			t.Errorf("Z1 = %.6f; want %.3f", got, wantZ1)
		}
	})

	t.Run("combined compressibility factor", func(t *testing.T) {
		got, err := CompressibilityFactor.Eval(Tr, Pr, w)
		if err != nil {
			t.Fatalf("Eval returned an unexpected error: %v", err)
		}

		// Looser than the table lookups: the published value combines
		// terms that were themselves rounded to three decimals.
		if math.Abs(got-wantZ) > 2e-3 {
			t.Errorf("Z = %.6f; want %.3f", got, wantZ)
		}
	})

	t.Run("molar volume", func(t *testing.T) {
		const (
			R      = 83.14 // bar·cm³/(mol·K)
			T      = 510.0
			P      = 25.0
			wantV  = 1480.7
			relTol = 2e-3
		)

		z, err := CompressibilityFactor.Eval(Tr, Pr, w)
		if err != nil {
			t.Fatalf("Eval returned an unexpected error: %v", err)
		}

		got := z * R * T / P

		if rel := math.Abs(got-wantV) / wantV; rel > relTol {
			t.Errorf("V = %.2f cm³/mol; want %.1f (%.3f%% apart)", got, wantV, 100*rel)
		}
	})
}

// TestTableStructure checks the invariants the interpolation relies on.
//
// Values is indexed as Values[TrIndex][PrIndex], so its outer length
// must match the temperature axis and every row must match the pressure
// axis. Both axes must be strictly increasing: the lookup uses a binary
// search, and the interpolation weights divide by the spacing between
// neighbouring points, which would be zero or negative otherwise.
//
// Each axis must also carry at least two points. findIndex reports -1
// for shorter axes and interpolate does not screen for it, so a
// degenerate table would index out of range. The tables are generated
// rather than supplied by callers, and the type is unexported, so the
// condition cannot be reached from outside the package — this pins the
// precondition rather than exercising a reachable path.
func TestTableStructure(t *testing.T) {
	for name, tbl := range tables() {
		t.Run(name, func(t *testing.T) {
			if len(tbl.Pr) < 2 {
				t.Fatalf("pressure axis has %d points; interpolation requires at least 2", len(tbl.Pr))
			}

			if len(tbl.Tr) < 2 {
				t.Fatalf("temperature axis has %d points; interpolation requires at least 2", len(tbl.Tr))
			}

			if len(tbl.Values) != len(tbl.Tr) {
				t.Fatalf(
					"Values has %d rows but the temperature axis has %d points",
					len(tbl.Values), len(tbl.Tr),
				)
			}

			for j, row := range tbl.Values {
				if len(row) != len(tbl.Pr) {
					t.Errorf(
						"row %d has %d entries but the pressure axis has %d points",
						j, len(row), len(tbl.Pr),
					)
				}
			}

			for i := 1; i < len(tbl.Pr); i++ {
				if tbl.Pr[i] <= tbl.Pr[i-1] {
					t.Errorf(
						"pressure axis is not strictly increasing at index %d: %g follows %g",
						i, tbl.Pr[i], tbl.Pr[i-1],
					)
				}
			}

			for j := 1; j < len(tbl.Tr); j++ {
				if tbl.Tr[j] <= tbl.Tr[j-1] {
					t.Errorf(
						"temperature axis is not strictly increasing at index %d: %g follows %g",
						j, tbl.Tr[j], tbl.Tr[j-1],
					)
				}
			}
		})
	}
}

// TestInterpolationAtGridPoints checks that the interpolation reproduces
// the tabulated values exactly at the nodes.
//
// Bilinear interpolation degenerates to a lookup when the query lands on
// a grid point, so any deviation indicates the weights or the index
// arithmetic are wrong. Evaluating every node of every table also
// confirms the two axes are not transposed, which a single sample near
// the middle of a square table could miss.
func TestInterpolationAtGridPoints(t *testing.T) {
	const tol = 1e-12

	for name, tbl := range tables() {
		t.Run(name, func(t *testing.T) {
			for j, tr := range tbl.Tr {
				for i, pr := range tbl.Pr {
					got, err := tbl.At(tr, pr)
					if err != nil {
						t.Fatalf("At(%g, %g) returned an unexpected error: %v", tr, pr, err)
					}

					if want := tbl.Values[j][i]; math.Abs(got-want) > tol {
						t.Errorf(
							"At(Tr=%g, Pr=%g) = %g; want the tabulated value %g",
							tr, pr, got, want,
						)
					}
				}
			}
		})
	}
}

// TestInterpolationAtCellCentres checks the interpolation against a
// closed-form property of the bilinear scheme: at the centre of a cell
// all four corners carry equal weight, so the result is their mean.
//
// This is independent of the tabulated data and fails if the weights are
// transposed or misnormalised — an error that grid points alone cannot
// reveal, since the weights collapse to zero and one there.
func TestInterpolationAtCellCentres(t *testing.T) {
	const tol = 1e-12

	for name, tbl := range tables() {
		t.Run(name, func(t *testing.T) {
			for j := 0; j < len(tbl.Tr)-1; j++ {
				for i := 0; i < len(tbl.Pr)-1; i++ {
					tr := (tbl.Tr[j] + tbl.Tr[j+1]) / 2
					pr := (tbl.Pr[i] + tbl.Pr[i+1]) / 2

					got, err := tbl.At(tr, pr)
					if err != nil {
						t.Fatalf("At(%g, %g) returned an unexpected error: %v", tr, pr, err)
					}

					want := (tbl.Values[j][i] +
						tbl.Values[j][i+1] +
						tbl.Values[j+1][i] +
						tbl.Values[j+1][i+1]) / 4

					if math.Abs(got-want) > tol {
						t.Errorf(
							"At the centre of cell (%d, %d): got %g; want the corner mean %g",
							j, i, got, want,
						)
					}
				}
			}
		})
	}
}

// TestInterpolationOutOfRange checks that conditions beyond the tabulated
// range are refused.
//
// The correlation carries no information outside its grid, so
// extrapolating would return a plausible number with no basis. Both ends
// of both axes are checked.
func TestInterpolationOutOfRange(t *testing.T) {
	for name, tbl := range tables() {
		t.Run(name, func(t *testing.T) {
			var (
				loPr = tbl.Pr[0]
				hiPr = tbl.Pr[len(tbl.Pr)-1]
				loTr = tbl.Tr[0]
				hiTr = tbl.Tr[len(tbl.Tr)-1]
			)

			testCases := []struct {
				name   string
				tr, pr float64
			}{
				{"pressure below range", loTr, loPr / 2},
				{"pressure above range", loTr, hiPr * 2},
				{"temperature below range", loTr / 2, loPr},
				{"temperature above range", hiTr * 2, loPr},
			}

			for _, tc := range testCases {
				t.Run(tc.name, func(t *testing.T) {
					if _, err := tbl.At(tc.tr, tc.pr); err == nil {
						t.Errorf("At(%g, %g): expected an error; got nil", tc.tr, tc.pr)
					}
				})
			}
		})
	}
}

// TestPropertyCombination checks that each property folds its two tables
// in the way its definition requires.
//
// The compressibility factor and the residual properties are linear in
// the acentric factor,
//
//	M = M⁰ + ω M¹,
//
// while the fugacity coefficient combines as a power law,
//
//	φ = φ⁰ (φ¹)^ω,
//
// because the linear relation applies to its logarithm. Mixing the two
// yields a plausible number, so the distinction is checked explicitly.
func TestPropertyCombination(t *testing.T) {
	const (
		Tr  = 1.200
		Pr  = 0.659
		w   = 0.200
		tol = 1e-12
	)

	additiveCases := []struct {
		name            string
		property        Property
		base, departure *table
	}{
		{"compressibility factor", CompressibilityFactor, Z0Table, Z1Table},
		{"residual enthalpy", ResidualEnthalpy, H0Table, H1Table},
		{"residual entropy", ResidualEntropy, S0Table, S1Table},
	}

	for _, tc := range additiveCases {
		t.Run(tc.name, func(t *testing.T) {
			m0, err := tc.base.At(Tr, Pr)
			if err != nil {
				t.Fatalf("base lookup failed: %v", err)
			}

			m1, err := tc.departure.At(Tr, Pr)
			if err != nil {
				t.Fatalf("departure lookup failed: %v", err)
			}

			got, err := tc.property.Eval(Tr, Pr, w)
			if err != nil {
				t.Fatalf("Eval returned an unexpected error: %v", err)
			}

			if want := m0 + w*m1; math.Abs(got-want) > tol {
				t.Errorf("Eval = %.12f; want m0 + w*m1 = %.12f", got, want)
			}
		})
	}

	t.Run("fugacity coefficient", func(t *testing.T) {
		p0, err := PHI0Table.At(Tr, Pr)
		if err != nil {
			t.Fatalf("base lookup failed: %v", err)
		}

		p1, err := PHI1Table.At(Tr, Pr)
		if err != nil {
			t.Fatalf("departure lookup failed: %v", err)
		}

		got, err := FugacityCoefficient.Eval(Tr, Pr, w)
		if err != nil {
			t.Fatalf("Eval returned an unexpected error: %v", err)
		}

		if want := p0 * math.Pow(p1, w); math.Abs(got-want) > tol {
			t.Errorf("Eval = %.12f; want p0 * p1^w = %.12f", got, want)
		}

		// The power law must differ from the additive form here, so the
		// check above is meaningful rather than coincidental.
		if additive := p0 + w*p1; math.Abs(got-additive) < 1e-6 {
			t.Errorf(
				"the power law and additive forms agree at this state (%.12f); choose another to make the test discriminating",
				got,
			)
		}
	})
}

// TestPropertyAtZeroAcentricFactor checks that a simple fluid reduces to
// the base term alone, since the departure term is weighted by the
// acentric factor.
func TestPropertyAtZeroAcentricFactor(t *testing.T) {
	const (
		Tr  = 1.200
		Pr  = 0.659
		tol = 1e-12
	)

	testCases := []struct {
		name     string
		property Property
		base     *table
	}{
		{"compressibility factor", CompressibilityFactor, Z0Table},
		{"residual enthalpy", ResidualEnthalpy, H0Table},
		{"residual entropy", ResidualEntropy, S0Table},
		{"fugacity coefficient", FugacityCoefficient, PHI0Table},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			want, err := tc.base.At(Tr, Pr)
			if err != nil {
				t.Fatalf("base lookup failed: %v", err)
			}

			got, err := tc.property.Eval(Tr, Pr, 0)
			if err != nil {
				t.Fatalf("Eval returned an unexpected error: %v", err)
			}

			if math.Abs(got-want) > tol {
				t.Errorf("Eval with omega = 0 gives %.12f; want the base term %.12f", got, want)
			}
		})
	}
}

// TestUninitialisedProperty checks that the zero value is refused.
//
// Property is a struct holding two table pointers and a combining
// function, so a zero value carries no tables at all. Without the guard
// it would dereference a nil pointer.
func TestUninitialisedProperty(t *testing.T) {
	var property Property

	if _, err := property.Eval(1.2, 0.659, 0.2); err == nil {
		t.Error("expected an error; got nil")
	}
}

// TestIdealGasLimit checks the physical behaviour of the base
// compressibility table: at the lowest tabulated pressure a gas is
// nearly ideal, so Z⁰ approaches unity.
//
// Only supercritical temperatures are checked, since below the critical
// temperature the low-pressure entries describe a condensed phase.
func TestIdealGasLimit(t *testing.T) {
	const tol = 0.01

	lowest := Z0Table.Pr[0]

	for _, tr := range Z0Table.Tr {
		if tr < 1 {
			continue
		}

		got, err := Z0Table.At(tr, lowest)
		if err != nil {
			t.Fatalf("At(%g, %g) returned an unexpected error: %v", tr, lowest, err)
		}

		if math.Abs(got-1) > tol {
			t.Errorf(
				"at Tr = %g and Pr = %g: Z0 = %.4f; want approximately 1",
				tr, lowest, got,
			)
		}
	}
}
