package cp_test

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cp"
)

// TestRangeErrorAccompaniesResult checks the contract for temperatures
// outside a correlation's fitted range: the extrapolated value is
// returned together with a *RangeError rather than discarded.
//
// The correlations are polynomial fits that remain defined outside their
// range, and a caller may legitimately accept a small extrapolation —
// Example 6.5 of Smith, Van Ness & Abbott does exactly that, evaluating
// carbon dioxide 5 K below the fitted lower bound. Returning zero would
// force the caller to reimplement the integral to get a number the
// package already computed.
func TestRangeErrorAccompaniesResult(t *testing.T) {
	gas := cp.CarbonDioxideGas

	// 293.15 K lies just below the fitted lower bound of 298.15 K.
	state1 := zfactor.Args{T: 343.15, P: 150e5, R: zfactor.RSI}
	state2 := zfactor.Args{T: 293.15, P: 15e5, R: zfactor.RSI}

	testCases := []struct {
		name string
		fn   func(a, b zfactor.Args) (float64, error)
	}{
		{"enthalpy", gas.IdealGasEnthalpyChange},
		{"entropy", gas.IdealGasEntropyChange},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got, err := tc.fn(state1, state2)

			if err == nil {
				t.Fatal("expected a range error; got nil")
			}

			var rangeErr *cp.RangeError
			if !errors.As(err, &rangeErr) {
				t.Fatalf("error is not a *cp.RangeError: %v", err)
			}

			if rangeErr.T != 293.15 {
				t.Errorf("RangeError reports T = %g; want the offending temperature 293.15", rangeErr.T)
			}

			if rangeErr.TMin != gas.TMin || rangeErr.TMax != gas.TMax {
				t.Errorf(
					"RangeError reports the range [%g, %g]; want [%g, %g]",
					rangeErr.TMin, rangeErr.TMax, gas.TMin, gas.TMax,
				)
			}

			if got == 0 {
				t.Error("the extrapolated value was discarded; want it returned alongside the error")
			}
		})
	}
}

// TestExample6_5Step2 reproduces step 2 of Example 6.5 of Smith, Van Ness
// & Abbott: the ideal-gas enthalpy and entropy changes for carbon
// dioxide between 70°C and 20°C, with the pressure falling from 150 to
// 15 bar.
//
// Only the entropy change is compared with the published value. The
// example prints the first heat-capacity coefficient as 5.547, but that
// is inconsistent with the tabulated Cp298/R of 4.467, which requires
// 5.457 — the value this package carries and the one its own entropy
// result agrees with. The published enthalpy change is not reproducible
// from either coefficient, so asserting it would pin a number neither
// the data nor the example supports.
func TestExample6_5Step2(t *testing.T) {
	const (
		want   = 13.08 // J/(mol*K)
		relTol = 1e-3
	)

	gas := cp.CarbonDioxideGas

	state1 := zfactor.Args{T: 343.15, P: 150e5, R: zfactor.RSI}
	state2 := zfactor.Args{T: 293.15, P: 15e5, R: zfactor.RSI}

	got, err := gas.IdealGasEntropyChange(state1, state2)

	// A range error is expected here and does not invalidate the result.
	var rangeErr *cp.RangeError
	if err != nil && !errors.As(err, &rangeErr) {
		t.Fatalf("IdealGasEntropyChange returned an unexpected error: %v", err)
	}

	if rel := math.Abs(got-want) / math.Abs(want); rel > relTol {
		t.Errorf("entropy change = %.4f J/(mol*K); want %.2f (%.3f%% apart)", got, want, 100*rel)
	}
}

// TestInvalidInputReturnsZero checks that genuinely invalid input is
// still refused outright.
//
// Unlike a range violation, these leave no meaningful value to return: a
// non-positive temperature or pressure puts the integrals outside their
// domain, and conflicting gas constants leave no basis for scaling the
// result.
func TestInvalidInputReturnsZero(t *testing.T) {
	gas := cp.MethaneGas
	valid := zfactor.Args{T: 300, P: 1e5, R: zfactor.RSI}

	testCases := []struct {
		name  string
		other zfactor.Args
	}{
		{"non-positive temperature", zfactor.Args{T: 0, P: 1e5, R: zfactor.RSI}},
		{"negative temperature", zfactor.Args{T: -10, P: 1e5, R: zfactor.RSI}},
		{"non-positive pressure", zfactor.Args{T: 300, P: 0, R: zfactor.RSI}},
		{"conflicting gas constants", zfactor.Args{T: 300, P: 1e5, R: 83.14}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			for _, fn := range []struct {
				name string
				call func(a, b zfactor.Args) (float64, error)
			}{
				{"enthalpy", gas.IdealGasEnthalpyChange},
				{"entropy", gas.IdealGasEntropyChange},
			} {
				got, err := fn.call(valid, tc.other)

				if err == nil {
					t.Errorf("%s: expected an error; got nil", fn.name)
					continue
				}

				var rangeErr *cp.RangeError
				if errors.As(err, &rangeErr) {
					t.Errorf("%s: invalid input reported as a range error", fn.name)
				}

				if got != 0 {
					t.Errorf("%s: got %g; want 0 for invalid input", fn.name, got)
				}
			}
		})
	}
}

// TestWithinRangeReportsNoError checks that a calculation wholly inside
// the fitted interval carries no error at all.
func TestWithinRangeReportsNoError(t *testing.T) {
	gas := cp.CarbonDioxideGas

	state1 := zfactor.Args{T: 400, P: 1e5, R: zfactor.RSI}
	state2 := zfactor.Args{T: 600, P: 1e5, R: zfactor.RSI}

	if _, err := gas.IdealGasEnthalpyChange(state1, state2); err != nil {
		t.Errorf("IdealGasEnthalpyChange returned an error inside the fitted range: %v", err)
	}

	if _, err := gas.IdealGasEntropyChange(state1, state2); err != nil {
		t.Errorf("IdealGasEntropyChange returned an error inside the fitted range: %v", err)
	}
}

// TestEnthalpyIntegralAgainstQuadrature checks the closed-form integral
// against numerical quadrature of the heat capacity it integrates.
//
// The antiderivative is written out by hand in the implementation, so
// nothing otherwise ties it to the Cp/R polynomial. Simpson's rule on a
// smooth quartic-free integrand converges quickly, so close agreement is
// expected.
func TestEnthalpyIntegralAgainstQuadrature(t *testing.T) {
	const (
		T1     = 400.0
		T2     = 900.0
		panels = 2000
		relTol = 1e-9
	)

	gas := cp.CarbonDioxideGas

	got, err := gas.IdealGasEnthalpyChange(
		zfactor.Args{T: T1, P: 1e5, R: zfactor.RSI},
		zfactor.Args{T: T2, P: 1e5, R: zfactor.RSI},
	)
	if err != nil {
		t.Fatalf("IdealGasEnthalpyChange returned an unexpected error: %v", err)
	}

	cpOverR := func(T float64) float64 {
		return gas.A + gas.B*T + gas.C*T*T + gas.D/(T*T)
	}

	want := zfactor.RSI * simpson(cpOverR, T1, T2, panels)

	if rel := math.Abs(got-want) / math.Abs(want); rel > relTol {
		t.Errorf("enthalpy change = %.9f J/mol; quadrature gives %.9f", got, want)
	}
}

// TestEntropyIntegralAgainstQuadrature applies the same check to the
// entropy integral, whose integrand is Cp/(R·T).
//
// The pressures are held equal so that only the temperature integral
// contributes.
func TestEntropyIntegralAgainstQuadrature(t *testing.T) {
	const (
		T1     = 400.0
		T2     = 900.0
		panels = 2000
		relTol = 1e-9
	)

	gas := cp.CarbonDioxideGas

	got, err := gas.IdealGasEntropyChange(
		zfactor.Args{T: T1, P: 1e5, R: zfactor.RSI},
		zfactor.Args{T: T2, P: 1e5, R: zfactor.RSI},
	)
	if err != nil {
		t.Fatalf("IdealGasEntropyChange returned an unexpected error: %v", err)
	}

	integrand := func(T float64) float64 {
		return (gas.A + gas.B*T + gas.C*T*T + gas.D/(T*T)) / T
	}

	want := zfactor.RSI * simpson(integrand, T1, T2, panels)

	if rel := math.Abs(got-want) / math.Abs(want); rel > relTol {
		t.Errorf("entropy change = %.9f J/(mol*K); quadrature gives %.9f", got, want)
	}
}

// TestPressureTermIsIsothermal checks that at constant temperature the
// entropy change reduces to the ideal-gas pressure term, -R ln(P2/P1).
func TestPressureTermIsIsothermal(t *testing.T) {
	const (
		T   = 500.0
		P1  = 1e5
		P2  = 5e5
		tol = 1e-12
	)

	gas := cp.CarbonDioxideGas

	got, err := gas.IdealGasEntropyChange(
		zfactor.Args{T: T, P: P1, R: zfactor.RSI},
		zfactor.Args{T: T, P: P2, R: zfactor.RSI},
	)
	if err != nil {
		t.Fatalf("IdealGasEntropyChange returned an unexpected error: %v", err)
	}

	if want := -zfactor.RSI * math.Log(P2/P1); math.Abs(got-want) > tol {
		t.Errorf("entropy change = %.12f; want -R*ln(P2/P1) = %.12f", got, want)
	}
}

// simpson integrates f over [a, b] using the composite Simpson rule with
// the given even number of panels.
func simpson(f func(float64) float64, a, b float64, panels int) float64 {
	if panels%2 != 0 {
		panels++
	}

	h := (b - a) / float64(panels)
	sum := f(a) + f(b)

	for i := 1; i < panels; i++ {
		x := a + h*float64(i)
		if i%2 == 0 {
			sum += 2 * f(x)
		} else {
			sum += 4 * f(x)
		}
	}

	return sum * h / 3
}
