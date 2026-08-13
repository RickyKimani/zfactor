package virial_test

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/virial"
)

// R is the gas constant in bar·cm³/(mol·K); pressures are then in bar,
// temperatures in kelvin, B in cm³/mol and C in cm⁶/mol².
const R = 83.14

// isopropanol returns the state used throughout: isopropanol vapor at
// 200 °C, with second and third virial coefficients from the worked
// example in Smith, Van Ness & Abbott.
func isopropanol(pressure float64) zfactor.Args {
	return zfactor.Args{
		T: 473.15,
		P: pressure,
		R: R,
		B: -388.0,
		C: -26000.0,
	}
}

// TestTwoTermVolumeMatchesCompressibility checks the two functions of the
// two-term equation against each other.
//
// Both express the same truncation, one as a volume and one as a
// compressibility factor, so Z must equal PV/(RT) exactly:
//
//	V = RT/P + B  =>  PV/(RT) = 1 + BP/(RT).
//
// They are written separately in the package, so nothing but this
// comparison ties them together.
func TestTwoTermVolumeMatchesCompressibility(t *testing.T) {
	const tol = 1e-12

	for _, pressure := range []float64{0.1, 1, 5, 10, 15} {
		args := isopropanol(pressure)

		volume, err := virial.SolveForVolumeTwoTerm(args)
		if err != nil {
			t.Fatalf("SolveForVolumeTwoTerm returned an unexpected error: %v", err)
		}

		z, err := virial.CompressibilityTwoTerm(args)
		if err != nil {
			t.Fatalf("CompressibilityTwoTerm returned an unexpected error: %v", err)
		}

		if want := pressure * volume / (R * args.T); math.Abs(z-want) > tol {
			t.Errorf(
				"at %g bar: Z = %.12f but PV/(RT) = %.12f",
				pressure, z, want,
			)
		}
	}
}

// TestThreeTermRootsSatisfyTheEquation checks the roots returned by the
// three-term solver against the equation they are roots of.
//
// The solver rearranges Z = 1 + B/V + C/V^2 into a cubic and hands it to
// the root finder, so a returned volume must reproduce the same
// compressibility factor by both routes: from the equation of state and
// from PV/(RT). The check therefore exercises the rearrangement, the
// root finder and CompressibilityThreeTerm together.
//
// Only real roots are examined; a molar volume must be real to describe
// a state. The identity holds for every root the cubic has, including
// the small spurious one described in TestThreeTermHighPressureLimit.
func TestThreeTermRootsSatisfyTheEquation(t *testing.T) {
	const tol = 1e-9

	for _, pressure := range []float64{1, 5, 10, 20} {
		args := isopropanol(pressure)

		roots, err := virial.SolveForVolumeThreeTerm(args)
		if err != nil {
			t.Fatalf("SolveForVolumeThreeTerm returned an unexpected error: %v", err)
		}

		var realRoots int

		for _, root := range roots {
			if math.Abs(imag(root)) > 1e-9 {
				continue
			}

			volume := real(root)
			if volume <= 0 {
				continue
			}

			realRoots++

			z, err := virial.CompressibilityThreeTerm(volume, args)
			if err != nil {
				t.Fatalf("CompressibilityThreeTerm returned an unexpected error: %v", err)
			}

			want := pressure * volume / (R * args.T)

			if rel := math.Abs(z-want) / math.Max(math.Abs(want), 1); rel > tol {
				t.Errorf(
					"at %g bar: the root V = %.6f gives Z = %.10f from the series but PV/(RT) = %.10f",
					pressure, volume, z, want,
				)
			}
		}

		if realRoots == 0 {
			t.Errorf("at %g bar: no positive real root was returned", pressure)
		}
	}
}

// TestTwoTermApproachesIdealGas checks that both two-term functions
// return ideal-gas behaviour as the pressure falls.
//
// The correction to Z is proportional to pressure, so it vanishes in
// that limit and the molar volume tends to RT/P.
func TestTwoTermApproachesIdealGas(t *testing.T) {
	args := isopropanol(1e-6)

	z, err := virial.CompressibilityTwoTerm(args)
	if err != nil {
		t.Fatalf("CompressibilityTwoTerm returned an unexpected error: %v", err)
	}

	if math.Abs(z-1) > 1e-6 {
		t.Errorf("Z = %.12f at %g bar; want approximately 1", z, args.P)
	}

	volume, err := virial.SolveForVolumeTwoTerm(args)
	if err != nil {
		t.Fatalf("SolveForVolumeTwoTerm returned an unexpected error: %v", err)
	}

	ideal := R * args.T / args.P

	if rel := math.Abs(volume-ideal) / ideal; rel > 1e-6 {
		t.Errorf("V = %.4f cm³/mol; want approximately the ideal-gas volume %.4f", volume, ideal)
	}
}

// TestNegativeSecondCoefficientLowersCompressibility checks the sign of
// the correction.
//
// A negative second virial coefficient means attraction between pairs of
// molecules predominates, drawing the gas into a smaller volume than an
// ideal gas would occupy, so Z falls below one and the volume below
// RT/P. Both truncations must agree on that direction.
func TestNegativeSecondCoefficientLowersCompressibility(t *testing.T) {
	args := isopropanol(10)

	z, err := virial.CompressibilityTwoTerm(args)
	if err != nil {
		t.Fatalf("CompressibilityTwoTerm returned an unexpected error: %v", err)
	}

	if z >= 1 {
		t.Errorf("Z = %.6f with a negative B; want a value below 1", z)
	}

	volume, err := virial.SolveForVolumeTwoTerm(args)
	if err != nil {
		t.Fatalf("SolveForVolumeTwoTerm returned an unexpected error: %v", err)
	}

	if ideal := R * args.T / args.P; volume >= ideal {
		t.Errorf("V = %.4f cm³/mol; want less than the ideal-gas volume %.4f", volume, ideal)
	}
}

// TestTwoTermRejectsHighPressure checks the validity limit.
//
// Truncating the pressure series after the second coefficient is only
// justified at low density. Rather than return a number the truncation
// does not support, the two-term functions refuse above 15 bar; the
// three-term form, which retains one more coefficient, has no such
// limit.
func TestTwoTermRejectsHighPressure(t *testing.T) {
	t.Run("at the limit", func(t *testing.T) {
		if _, err := virial.CompressibilityTwoTerm(isopropanol(15)); err != nil {
			t.Errorf("15 bar should be accepted; got %v", err)
		}
	})

	t.Run("above the limit", func(t *testing.T) {
		if _, err := virial.CompressibilityTwoTerm(isopropanol(15.001)); err == nil {
			t.Error("expected an error above 15 bar; got nil")
		}

		if _, err := virial.SolveForVolumeTwoTerm(isopropanol(50)); err == nil {
			t.Error("expected an error above 15 bar; got nil")
		}
	})

	t.Run("three-term has no such limit", func(t *testing.T) {
		if _, err := virial.SolveForVolumeThreeTerm(isopropanol(50)); err != nil {
			t.Errorf("the three-term equation should accept 50 bar; got %v", err)
		}
	})
}

// TestThreeTermHighPressureLimit records where the truncation stops
// describing a state, and what the solver does there.
//
// Two positive roots exist at low pressure: the larger is the physical
// vapor volume, near RT/P, and the smaller is an artefact of truncating
// the series, at a density where the neglected terms would dominate. As
// the pressure rises the two approach one another and, for isopropanol
// with these coefficients, merge and leave the real axis between 20 and
// 25 bar. Above that the equation has no real volume at all.
//
// The solver reports this by returning complex roots rather than an
// error, matching the cubic equation-of-state solver: selecting the
// physically meaningful root is the caller's part. A caller that wants
// the vapor volume must therefore take the largest real root and be
// prepared for there to be none.
func TestThreeTermHighPressureLimit(t *testing.T) {
	positiveRealRoots := func(pressure float64) []float64 {
		t.Helper()

		roots, err := virial.SolveForVolumeThreeTerm(isopropanol(pressure))
		if err != nil {
			t.Fatalf("SolveForVolumeThreeTerm returned an unexpected error: %v", err)
		}

		var real_ []float64
		for _, root := range roots {
			if math.Abs(imag(root)) < 1e-9 && real(root) > 0 {
				real_ = append(real_, real(root))
			}
		}

		return real_
	}

	t.Run("two positive roots below the limit", func(t *testing.T) {
		for _, pressure := range []float64{1, 10, 20} {
			got := positiveRealRoots(pressure)

			if len(got) != 2 {
				t.Errorf("at %g bar there are %d positive real roots; want 2", pressure, len(got))
				continue
			}

			// The larger root is the vapor volume and must sit near the
			// ideal-gas value; the smaller is the spurious one.
			ideal := R * 473.15 / pressure
			largest := math.Max(got[0], got[1])

			if largest > ideal {
				t.Errorf(
					"at %g bar the vapor volume %.2f exceeds the ideal-gas volume %.2f, though B is negative",
					pressure, largest, ideal,
				)
			}
		}
	})

	t.Run("no real root above the limit", func(t *testing.T) {
		for _, pressure := range []float64{25, 40, 50} {
			if got := positiveRealRoots(pressure); len(got) != 0 {
				t.Errorf("at %g bar the roots %v are real; the truncation is expected to have failed", pressure, got)
			}
		}
	})
}

// TestInvalidInput checks the guards on every entry point.
//
// A missing coefficient is rejected rather than silently treated as
// zero, since a zero coefficient would quietly reduce the equation to a
// lower truncation, or in the three-term case to one the solver cannot
// form a cubic from.
func TestInvalidInput(t *testing.T) {
	valid := isopropanol(10)

	with := func(mutate func(*zfactor.Args)) zfactor.Args {
		args := valid
		mutate(&args)
		return args
	}

	twoTermCases := []struct {
		name string
		args zfactor.Args
	}{
		{"zero pressure", with(func(a *zfactor.Args) { a.P = 0 })},
		{"negative pressure", with(func(a *zfactor.Args) { a.P = -1 })},
		{"zero temperature", with(func(a *zfactor.Args) { a.T = 0 })},
		{"negative temperature", with(func(a *zfactor.Args) { a.T = -10 })},
		{"zero gas constant", with(func(a *zfactor.Args) { a.R = 0 })},
		{"missing second coefficient", with(func(a *zfactor.Args) { a.B = 0 })},
	}

	for _, tc := range twoTermCases {
		t.Run("two-term/"+tc.name, func(t *testing.T) {
			if _, err := virial.SolveForVolumeTwoTerm(tc.args); err == nil {
				t.Error("SolveForVolumeTwoTerm: expected an error; got nil")
			}

			if _, err := virial.CompressibilityTwoTerm(tc.args); err == nil {
				t.Error("CompressibilityTwoTerm: expected an error; got nil")
			}
		})
	}

	threeTermCases := []struct {
		name string
		args zfactor.Args
	}{
		{"zero pressure", with(func(a *zfactor.Args) { a.P = 0 })},
		{"zero temperature", with(func(a *zfactor.Args) { a.T = 0 })},
		{"zero gas constant", with(func(a *zfactor.Args) { a.R = 0 })},
		{"missing second coefficient", with(func(a *zfactor.Args) { a.B = 0 })},
		{"missing third coefficient", with(func(a *zfactor.Args) { a.C = 0 })},
	}

	for _, tc := range threeTermCases {
		t.Run("three-term/"+tc.name, func(t *testing.T) {
			if _, err := virial.SolveForVolumeThreeTerm(tc.args); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}

	compressibilityCases := []struct {
		name   string
		volume float64
		args   zfactor.Args
	}{
		{"zero volume", 0, valid},
		{"negative volume", -100, valid},
		{"missing second coefficient", 3000, with(func(a *zfactor.Args) { a.B = 0 })},
		{"missing third coefficient", 3000, with(func(a *zfactor.Args) { a.C = 0 })},
	}

	for _, tc := range compressibilityCases {
		t.Run("three-term compressibility/"+tc.name, func(t *testing.T) {
			if _, err := virial.CompressibilityThreeTerm(tc.volume, tc.args); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestTruncationsAgreeAtLowDensity checks that the two truncations
// converge on each other where both are valid.
//
// The third coefficient enters through C/V^2, which becomes negligible
// against B/V as the gas thins, so the two- and three-term
// compressibility factors must approach one another at low pressure and
// separate as the density rises.
func TestTruncationsAgreeAtLowDensity(t *testing.T) {
	largestRealRoot := func(args zfactor.Args) float64 {
		t.Helper()

		roots, err := virial.SolveForVolumeThreeTerm(args)
		if err != nil {
			t.Fatalf("SolveForVolumeThreeTerm returned an unexpected error: %v", err)
		}

		largest := math.Inf(-1)
		for _, root := range roots {
			if math.Abs(imag(root)) < 1e-9 && real(root) > largest {
				largest = real(root)
			}
		}

		return largest
	}

	gap := func(pressure float64) float64 {
		args := isopropanol(pressure)

		twoTerm, err := virial.CompressibilityTwoTerm(args)
		if err != nil {
			t.Fatalf("CompressibilityTwoTerm returned an unexpected error: %v", err)
		}

		threeTerm, err := virial.CompressibilityThreeTerm(largestRealRoot(args), args)
		if err != nil {
			t.Fatalf("CompressibilityThreeTerm returned an unexpected error: %v", err)
		}

		return math.Abs(twoTerm - threeTerm)
	}

	dilute := gap(0.5)
	dense := gap(10)

	if dilute > 1e-3 {
		t.Errorf("the truncations differ by %.6f at 0.5 bar; want close agreement", dilute)
	}

	if dense <= dilute {
		t.Errorf(
			"the truncations differ by %.6f at 10 bar and %.6f at 0.5 bar; the gap should widen with density",
			dense, dilute,
		)
	}
}
