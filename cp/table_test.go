package cp_test

import (
	"math"
	"strings"
	"testing"

	"github.com/rickykimani/zfactor/cp"
)

// cp298Exceptions lists correlations whose tabulated Cp298/R disagrees
// with the constants published alongside it.
//
// Both are reproduced faithfully from Appendix C; the inconsistency is
// in the source. In each case the constants appear to be the reliable
// pair and the convenience column the outlier: evaluating them at
// 298.15 K lands nearer the accepted heat capacity than the tabulated
// Cp298/R does.
// The key is the name together with the formula, which distinguishes
// entries appearing in more than one table. 1,3-Butadiene is listed as
// both a gas and a liquid, and only the gas is inconsistent: the
// liquid's constants reproduce its tabulated Cp298 exactly. Keying on
// the name alone would excuse both.
var cp298Exceptions = map[string]string{
	"1,3-Butadiene/C4H6": "Table C.1 gives Cp298/R = 10.720 where its own constants give 9.931",
	"S (rhombic)/":       "Table C.2 gives Cp298/R = 3.748 where its own constants give 2.718",
}

// phaseKey identifies an entry across the three tables. Liquids carry no
// formula in the source, so theirs is empty.
func phaseKey(entry *cp.HeatCapacity) string {
	return entry.Name + "/" + entry.Formula
}

// TestCp298Consistency checks each tabulated Cp298/R against the
// correlation published with it,
//
//	Cp/R = A + B*T + C*T^2 + D*T^-2   at T = 298.15 K.
//
// Cp298 is a convenience value rather than an independent measurement,
// so the two must agree. The check exists to guard the data pipeline
// rather than the thermodynamics: the parser reads fixed-width columns
// out of the textbook's PDF, and a misread shifts every later value one
// place along the row.
//
// That is not hypothetical. Chemical subscripts are rendered as separate
// text runs, so NH3 extracted as "NH" followed by "3"; the parser took
// the stray digit for the first data column and shifted the rest, which
// silently corrupted eight entries and dropped their final coefficient.
// This check identifies every one of them.
func TestCp298Consistency(t *testing.T) {
	const (
		T   = 298.15
		tol = 5e-3
	)

	if len(cp.All) == 0 {
		t.Fatal("cp.All is empty")
	}

	for _, entry := range cp.All {
		t.Run(phaseKey(entry), func(t *testing.T) {
			if reason, skip := cp298Exceptions[phaseKey(entry)]; skip {
				t.Skipf("known source inconsistency: %s", reason)
			}

			if entry.Cp298 == 0 {
				t.Skip("no tabulated Cp298")
			}

			got := entry.A + entry.B*T + entry.C*T*T + entry.D/(T*T)

			if rel := math.Abs(got-entry.Cp298) / entry.Cp298; rel > tol {
				t.Errorf(
					"tabulated Cp298/R = %.4f but the constants give %.4f (%.2f%% apart)",
					entry.Cp298, got, 100*rel,
				)
			}
		})
	}
}

// TestTemperatureRangesAreSensible checks that every correlation covers
// a usable interval.
//
// The corrupted entries carried a TMax of 2 or 3 K, inherited from a
// formula subscript. That rejected every temperature a caller could
// supply, which is precisely why the corruption went unnoticed: the
// package never returned a wrong number, it returned nothing at all.
func TestTemperatureRangesAreSensible(t *testing.T) {
	for _, entry := range cp.All {
		t.Run(phaseKey(entry), func(t *testing.T) {
			if entry.TMax <= entry.TMin {
				t.Errorf(
					"TMax = %g does not exceed TMin = %g",
					entry.TMax, entry.TMin,
				)
			}

			// Every correlation in Appendix C is fitted from 298.15 K
			// upward, so a range that excludes room temperature means
			// the columns have shifted.
			if entry.TMax < 298.15 {
				t.Errorf(
					"TMax = %g K lies below the fitted lower bound of 298.15 K",
					entry.TMax,
				)
			}
		})
	}
}

// TestFormulasKeepTheirSubscripts checks that gas formulas were not
// truncated when the table was parsed.
//
// A formula ending in a letter where the species name implies a
// subscript — "NH" for ammonia — is the visible symptom of the column
// shift that TestCp298Consistency detects numerically.
func TestFormulasKeepTheirSubscripts(t *testing.T) {
	testCases := []struct {
		entry *cp.HeatCapacity
		want  string
	}{
		{cp.AmmoniaGas, "NH3"},
		{cp.NitrogenDioxideGas, "NO2"},
		{cp.MethaneGas, "CH4"},
		{cp.CarbonDioxideGas, "CO2"},
	}

	for _, tc := range testCases {
		t.Run(tc.want, func(t *testing.T) {
			if got := strings.TrimSpace(tc.entry.Formula); got != tc.want {
				t.Errorf("formula = %q; want %q", got, tc.want)
			}
		})
	}
}

// TestSolidNamesKeepTheirSubscripts checks the same for solids, whose
// formula doubles as the species name in Table C.2.
//
// A subscript that splits mid-formula leaves a space rather than
// truncating the name — NH4Cl became "NH 4Cl" — so this catches a class
// of parsing defect the numeric checks cannot: the values are correct,
// only the label is wrong.
func TestSolidNamesKeepTheirSubscripts(t *testing.T) {
	testCases := []struct {
		entry *cp.HeatCapacity
		want  string
	}{
		{cp.NH4ClSolid, "NH4Cl"},
		{cp.Fe2O3Solid, "Fe2O3"},
		{cp.Fe3O4Solid, "Fe3O4"},
		{cp.I2Solid, "I2"},
		{cp.CaCO3Solid, "CaCO3"},
		{cp.CaOH2Solid, "Ca(OH)2"},
		{cp.SiO2QuartzSolid, "SiO2 (quartz)"},
	}

	for _, tc := range testCases {
		t.Run(tc.want, func(t *testing.T) {
			if got := tc.entry.Name; got != tc.want {
				t.Errorf("name = %q; want %q", got, tc.want)
			}
		})
	}
}

// TestCorrectedEntries pins the values of the entries repaired after the
// subscript-shift bug, so a regression in the parser is caught here
// rather than in a calculation far downstream.
//
// The values are those of Appendix C, each verified against the Cp298
// identity.
func TestCorrectedEntries(t *testing.T) {
	testCases := []struct {
		name                string
		entry               *cp.HeatCapacity
		tMax, cp298         float64
		a, bScaled, dScaled float64
	}{
		{"Ammonia", cp.AmmoniaGas, 1800, 4.269, 3.578, 3.020, -0.186},
		{"Nitrogen dioxide", cp.NitrogenDioxideGas, 2000, 4.447, 4.982, 1.195, -0.792},
		{"CaCO3", cp.CaCO3Solid, 1200, 9.848, 12.572, 2.637, -3.120},
		{"CaC2", cp.CaC2Solid, 720, 7.508, 8.254, 1.429, -1.042},
		{"CaCl2", cp.CaCl2Solid, 1055, 8.762, 8.646, 1.530, -0.302},
		{"SiO2 (quartz)", cp.SiO2QuartzSolid, 847, 5.345, 4.871, 5.365, -1.001},
	}

	const tol = 1e-9

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if tc.entry.TMax != tc.tMax {
				t.Errorf("TMax = %g; want %g", tc.entry.TMax, tc.tMax)
			}

			if tc.entry.Cp298 != tc.cp298 {
				t.Errorf("Cp298 = %g; want %g", tc.entry.Cp298, tc.cp298)
			}

			if tc.entry.A != tc.a {
				t.Errorf("A = %g; want %g", tc.entry.A, tc.a)
			}

			if math.Abs(tc.entry.B-tc.bScaled*1e-3) > tol {
				t.Errorf("B = %g; want %g", tc.entry.B, tc.bScaled*1e-3)
			}

			if math.Abs(tc.entry.D-tc.dScaled*1e5) > tol*1e5 {
				t.Errorf("D = %g; want %g", tc.entry.D, tc.dScaled*1e5)
			}
		})
	}
}
