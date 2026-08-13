package cubic_test

import (
	"math"
	"math/cmplx"
	"math/rand"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/substance"
)

// R is the gas constant in bar·cm³/(mol·K), matching the units used by
// the worked examples: temperatures in kelvin, pressures in bar and
// molar volumes in cm³/mol.
const R = 10 * zfactor.RSI

// TestSolveForVolumeExample3_9 reproduces Example 3.9 of Smith, Van Ness
// & Abbott: the saturated-vapor and saturated-liquid molar volumes of
// n-butane at 350 K, where the vapor pressure is 9.4573 bar.
//
// The example is evaluated at saturation, so the cubic has three real
// roots for every equation of state. The smallest is the saturated
// liquid, the largest the saturated vapor, and the middle root lies on
// the mechanically unstable branch and has no physical meaning.
//
// The published comparison table gives both volumes for all four cubic
// equations, which makes this a single reference point covering the
// entire EOS family.
func TestSolveForVolumeExample3_9(t *testing.T) {
	const (
		T = 350.0  // K
		P = 9.4573 // bar, the vapor pressure at 350 K
		// Relative tolerance. The published volumes carry four
		// significant figures, and every equation of state reproduces
		// them to within a few hundredths of a percent, so 0.1% is
		// comfortably tight enough to detect an incorrect parameter.
		tol = 1e-3
	)

	testCases := []struct {
		name      string
		eos       cubic.EOSType
		wantLiq   float64
		wantVapor float64
	}{
		{"vdW", &cubic.VdW{}, 191.0, 2667},
		{"RK", &cubic.RK{}, 133.3, 2555},
		{"SRK", &cubic.SRK{}, 127.8, 2520},
		{"PR", &cubic.PR{}, 112.6, 2486},
	}

	nButane := substance.NButane
	args := zfactor.Args{T: T, P: P, R: R}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			cfg := nButane.CubicConfig(tc.eos, args)

			res, err := cubic.SolveForVolume(cfg)
			if err != nil {
				t.Fatalf("SolveForVolume returned an unexpected error: %v", err)
			}

			roots := res.Clean()
			if len(roots) != 3 {
				t.Fatalf(
					"expected three real roots at saturation; got %d: %v",
					len(roots), roots,
				)
			}

			liquid, middle, vapor := roots[0], roots[1], roots[2]

			if relErr(liquid, tc.wantLiq) > tol {
				t.Errorf(
					"saturated-liquid volume = %.2f cm³/mol; want %.1f (%.2f%% off)",
					liquid, tc.wantLiq, 100*relErr(liquid, tc.wantLiq),
				)
			}

			if relErr(vapor, tc.wantVapor) > tol {
				t.Errorf(
					"saturated-vapor volume = %.1f cm³/mol; want %.0f (%.2f%% off)",
					vapor, tc.wantVapor, 100*relErr(vapor, tc.wantVapor),
				)
			}

			// The unstable root must separate the two physical roots.
			if middle <= liquid || middle >= vapor {
				t.Errorf(
					"middle root %.2f must lie between the liquid (%.2f) and vapor (%.2f) roots",
					middle, liquid, vapor,
				)
			}

			// Every root must satisfy the equation of state it came from.
			for _, v := range roots {
				pres, err := cubic.Pressure(cfg, v)
				if err != nil {
					t.Fatalf("Pressure returned an unexpected error: %v", err)
				}
				if relErr(pres.P, P) > 1e-6 {
					t.Errorf(
						"root V = %.4f does not reproduce the specified pressure: got %.6f bar, want %.4f",
						v, pres.P, P,
					)
				}
			}
		})
	}
}

// TestEOSParameters checks the equation-of-state constants against their
// closed-form definitions.
//
// For the Redlich/Kwong family the constants follow from the critical
// isotherm conditions and are exact:
//
//	Ω = (2^(1/3) - 1)/3        = 0.08664
//	Ψ = 1/(9(2^(1/3) - 1))     = 0.42748
//
// The van der Waals constants are exact rationals. The Peng/Robinson
// constants are the published values, which are themselves rounded
// roots of the corresponding critical conditions and so are compared
// only to the precision at which they are quoted.
func TestEOSParameters(t *testing.T) {
	c := math.Cbrt(2) - 1

	testCases := []struct {
		name       string
		eos        cubic.EOSType
		omega, psi float64
		tol        float64
	}{
		{"vdW", &cubic.VdW{}, 1.0 / 8.0, 27.0 / 64.0, 1e-12},
		{"RK", &cubic.RK{}, c / 3, 1 / (9 * c), 1e-5},
		{"SRK", &cubic.SRK{}, c / 3, 1 / (9 * c), 1e-5},
		{"PR", &cubic.PR{}, 0.07780, 0.45724, 1e-12},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got := tc.eos.Params()

			if math.Abs(got.Omega-tc.omega) > tc.tol {
				t.Errorf("Omega = %.7f; want %.7f", got.Omega, tc.omega)
			}

			if math.Abs(got.Psi-tc.psi) > tc.tol {
				t.Errorf("Psi = %.7f; want %.7f", got.Psi, tc.psi)
			}
		})
	}
}

// TestAlphaAtCriticalTemperature checks that every alpha function
// reduces to unity at the critical temperature, where Tr = 1.
//
// The Soave-type expressions are built as [1 + m(1 - sqrt(Tr))]², so
// this holds for any acentric factor and is independent of the fitted
// coefficients.
func TestAlphaAtCriticalTemperature(t *testing.T) {
	testCases := []struct {
		name string
		eos  cubic.EOSType
	}{
		{"vdW", &cubic.VdW{}},
		{"RK", &cubic.RK{}},
		{"SRK", &cubic.SRK{}},
		{"PR", &cubic.PR{}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			for _, w := range []float64{0, 0.2, 0.35} {
				if got := tc.eos.Alpha(1, w); math.Abs(got-1) > 1e-12 {
					t.Errorf("Alpha(1, %.2f) = %.12f; want 1", w, got)
				}
			}
		})
	}
}

// TestSolveForVolumeSupercritical checks that a state well above the
// critical point yields exactly one real root, since no phase boundary
// exists there.
func TestSolveForVolumeSupercritical(t *testing.T) {
	co2 := substance.CarbonDioxide

	// Tr = 1.5, Pr = 2, comfortably supercritical.
	args := zfactor.Args{
		T: 1.5 * co2.Critical.Tc,
		P: 2.0 * co2.Critical.Pc,
		R: R,
	}

	for _, eos := range []cubic.EOSType{
		&cubic.VdW{}, &cubic.RK{}, &cubic.SRK{}, &cubic.PR{},
	} {
		cfg := co2.CubicConfig(eos, args)

		res, err := cubic.SolveForVolume(cfg)
		if err != nil {
			t.Fatalf("SolveForVolume returned an unexpected error: %v", err)
		}

		roots := res.Clean()
		if len(roots) != 1 {
			t.Errorf(
				"supercritical state should yield one real root; got %d: %v",
				len(roots), roots,
			)
		}
	}
}

// TestSolveForVolumeInvalidInput checks that invalid states are
// rejected rather than silently producing roots.
func TestSolveForVolumeInvalidInput(t *testing.T) {
	testCases := []struct {
		name string
		cfg  *cubic.EOSCfg
	}{
		{
			name: "non-positive temperature",
			cfg:  &cubic.EOSCfg{Type: &cubic.PR{}, T: 0, P: 10, Tc: 425.1, Pc: 37.96, R: R},
		},
		{
			name: "non-positive pressure",
			cfg:  &cubic.EOSCfg{Type: &cubic.PR{}, T: 350, P: 0, Tc: 425.1, Pc: 37.96, R: R},
		},
		{
			name: "invalid critical property",
			cfg:  &cubic.EOSCfg{Type: &cubic.PR{}, T: 350, P: 10, Tc: 0, Pc: 37.96, R: R},
		},
		{
			name: "invalid gas constant",
			cfg:  &cubic.EOSCfg{Type: &cubic.PR{}, T: 350, P: 10, Tc: 425.1, Pc: 37.96, R: 0},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := cubic.SolveForVolume(tc.cfg); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestSolveCubicResidual is a property test over the root solver
// underlying every equation of state.
//
// Rather than comparing against known roots, it substitutes each
// returned root back into the polynomial and requires the result to
// vanish. This holds for any cubic with any root structure, so it
// covers regions no worked example visits — in particular the
// single-real-root case, where an earlier implementation returned roots
// with residuals of order 10⁶.
//
// The residual is scaled by the magnitude of the individual terms so
// that the criterion remains meaningful for poorly scaled coefficients.
//
// The polynomials are monic, matching those SolveForVolume constructs:
// its leading coefficient is always 1. Cardano's method loses accuracy
// when the leading coefficient is small enough that the cubic is nearly
// a quadratic, since one root then diverges and the depressed-cubic
// shift cancels catastrophically. That regime is unreachable through
// the equation-of-state solvers and is excluded here.
//
// The threshold is set well inside the observed accuracy: across a wide
// sweep of equation-of-state states the worst relative error is of
// order 1e-5, while a genuinely incorrect root — such as those produced
// by the earlier cube-root branch — has a relative residual of order 1.
func TestSolveCubicResidual(t *testing.T) {
	const tol = 1e-4

	rng := rand.New(rand.NewSource(1))

	for i := 0; i < 5000; i++ {
		a := 1.0
		b := coefficient(rng)
		c := coefficient(rng)
		d := coefficient(rng)

		roots, err := zfactor.SolveCubic(a, b, c, d)
		if err != nil {
			t.Fatalf("SolveCubic(%g, %g, %g, %g) returned an unexpected error: %v",
				a, b, c, d, err)
		}

		for j, x := range roots {
			residual, scale := evaluateCubic(a, b, c, d, x)

			if cmplx.Abs(residual) > tol*scale {
				t.Fatalf(
					"SolveCubic(%g, %g, %g, %g): root %d = %v has residual %g (scale %g)",
					a, b, c, d, j, x, cmplx.Abs(residual), scale,
				)
			}
		}
	}
}

// TestSolveCubicConjugatePairs checks a structural property of cubics
// with real coefficients: the roots are either all real, or one real
// root together with a complex-conjugate pair. A set with no real root,
// or with unmatched imaginary parts, is impossible.
func TestSolveCubicConjugatePairs(t *testing.T) {
	rng := rand.New(rand.NewSource(2))

	for i := 0; i < 5000; i++ {
		a := 1.0
		b := coefficient(rng)
		c := coefficient(rng)
		d := coefficient(rng)

		roots, err := zfactor.SolveCubic(a, b, c, d)
		if err != nil {
			t.Fatalf("unexpected error: %v", err)
		}

		// Scale the "is real" threshold to the size of the roots.
		var largest float64
		for _, x := range roots {
			largest = math.Max(largest, cmplx.Abs(x))
		}
		tol := 1e-9 * math.Max(largest, 1)

		var real, complexRoots []complex128
		for _, x := range roots {
			if math.Abs(imag(x)) <= tol {
				real = append(real, x)
			} else {
				complexRoots = append(complexRoots, x)
			}
		}

		if len(real) == 0 {
			t.Fatalf(
				"SolveCubic(%g, %g, %g, %g) returned no real root: %v",
				a, b, c, d, roots,
			)
		}

		if len(complexRoots) == 1 {
			t.Fatalf(
				"SolveCubic(%g, %g, %g, %g) returned an unpaired complex root: %v",
				a, b, c, d, roots,
			)
		}

		if len(complexRoots) == 2 {
			sum := imag(complexRoots[0]) + imag(complexRoots[1])
			if math.Abs(sum) > tol {
				t.Fatalf(
					"SolveCubic(%g, %g, %g, %g): complex roots are not conjugates: %v",
					a, b, c, d, complexRoots,
				)
			}
		}
	}
}

// evaluateCubic returns the value of ax³ + bx² + cx + d at x together
// with the magnitude of its largest term, used to scale the residual.
func evaluateCubic(a, b, c, d float64, x complex128) (complex128, float64) {
	x2 := x * x
	x3 := x2 * x

	terms := []complex128{
		complex(a, 0) * x3,
		complex(b, 0) * x2,
		complex(c, 0) * x,
		complex(d, 0),
	}

	var value complex128
	scale := 1.0
	for _, term := range terms {
		value += term
		scale = math.Max(scale, cmplx.Abs(term))
	}

	return value, scale
}

// coefficient returns a random polynomial coefficient spanning several
// orders of magnitude, matching the wide dynamic range seen in
// equation-of-state coefficients.
func coefficient(rng *rand.Rand) float64 {
	return (rng.Float64()*2 - 1) * math.Pow(10, rng.Float64()*6-2)
}

// relErr returns the relative difference between got and want.
func relErr(got, want float64) float64 {
	return math.Abs(got-want) / math.Abs(want)
}
