package leekesler

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
)

// n-butane, used throughout: the critical constants and normal boiling
// point from Table B.1 of Smith, Van Ness & Abbott.
const (
	butaneTn = 272.70
	butaneTc = 425.10
	butanePc = 37.96
	butaneW  = 0.200
)

// TestVaporPressureAtNormalBoilingPoint checks the correlation against
// its own definition.
//
// The normal boiling point is where the saturation pressure equals one
// atmosphere, and VaporPressure estimates the acentric factor from that
// same point before evaluating. Supplying Tn as the temperature must
// therefore return one atmosphere, whatever the substance and however
// the correlation performs elsewhere.
func TestVaporPressureAtNormalBoilingPoint(t *testing.T) {
	const tol = 1e-9

	got, err := VaporPressure(butaneTn, butaneTn, butaneTc, butanePc)
	if err != nil {
		t.Fatalf("VaporPressure returned an unexpected error: %v", err)
	}

	if rel := math.Abs(got-zfactor.AtmBar) / zfactor.AtmBar; rel > tol {
		t.Errorf(
			"vapor pressure at the normal boiling point = %.9f bar; want %.9f bar (1 atm)",
			got, zfactor.AtmBar,
		)
	}
}

// TestVaporPressureAtCriticalTemperature checks that the correlation
// returns approximately the critical pressure at the critical
// temperature, where the reduced vapor pressure reaches unity.
//
// The agreement is approximate because the two reference-fluid
// expressions are empirical fits rather than exact relations, so the
// tolerance reflects the correlation's own accuracy.
func TestVaporPressureAtCriticalTemperature(t *testing.T) {
	const tol = 0.02

	got, err := VaporPressure(butaneTc, butaneTn, butaneTc, butanePc)
	if err != nil {
		t.Fatalf("VaporPressure returned an unexpected error: %v", err)
	}

	if rel := math.Abs(got-butanePc) / butanePc; rel > tol {
		t.Errorf(
			"vapor pressure at the critical temperature = %.4f bar; want approximately Pc = %.2f bar (%.2f%% apart)",
			got, butanePc, 100*rel,
		)
	}
}

// TestVaporPressureIncreasesWithTemperature checks that the saturation
// curve rises monotonically, as the Clausius-Clapeyron relation requires
// for any substance.
func TestVaporPressureIncreasesWithTemperature(t *testing.T) {
	var previous float64

	for _, tr := range []float64{0.6, 0.7, 0.8, 0.9, 0.95, 1.0} {
		got, err := VaporPressure(tr*butaneTc, butaneTn, butaneTc, butanePc)
		if err != nil {
			t.Fatalf("VaporPressure returned an unexpected error: %v", err)
		}

		if got <= previous {
			t.Errorf(
				"at Tr = %.2f the vapor pressure is %.6f bar; must exceed the previous value %.6f",
				tr, got, previous,
			)
		}

		previous = got
	}
}

// TestEstimateAcentricFactor checks the estimate against the tabulated
// acentric factor of n-butane.
//
// The correlation is fitted to normal fluids, of which n-butane is one,
// so the two should agree closely.
func TestEstimateAcentricFactor(t *testing.T) {
	const tol = 0.01

	got, err := EstimateAcentricFactor(butaneTn, butaneTc, butanePc)
	if err != nil {
		t.Fatalf("EstimateAcentricFactor returned an unexpected error: %v", err)
	}

	if math.Abs(got-butaneW) > tol {
		t.Errorf("acentric factor = %.4f; want approximately %.3f", got, butaneW)
	}
}

// TestVaporPressureInvalidInput checks that non-physical inputs are
// rejected rather than producing a number from a negative or zero
// reduced temperature.
func TestVaporPressureInvalidInput(t *testing.T) {
	testCases := []struct {
		name          string
		T, Tn, Tc, Pc float64
	}{
		{"non-positive temperature", 0, butaneTn, butaneTc, butanePc},
		{"non-positive critical temperature", 350, butaneTn, 0, butanePc},
		{"non-positive critical pressure", 350, butaneTn, butaneTc, 0},
		{"non-positive normal boiling point", 350, 0, butaneTc, butanePc},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := VaporPressure(tc.T, tc.Tn, tc.Tc, tc.Pc); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestEstimateAcentricFactorInvalidInput checks the same guards on the
// acentric-factor estimate.
func TestEstimateAcentricFactorInvalidInput(t *testing.T) {
	testCases := []struct {
		name       string
		Tn, Tc, Pc float64
	}{
		{"non-positive critical temperature", butaneTn, 0, butanePc},
		{"non-positive critical pressure", butaneTn, butaneTc, 0},
		{"non-positive normal boiling point", 0, butaneTc, butanePc},
		{"negative normal boiling point", -10, butaneTc, butanePc},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := EstimateAcentricFactor(tc.Tn, tc.Tc, tc.Pc); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}
