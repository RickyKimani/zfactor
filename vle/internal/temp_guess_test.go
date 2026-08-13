package internal

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/antoine"
)

// stubModel is an antoine.Model whose saturation temperature is fixed,
// so that the ordering logic can be exercised without depending on the
// tabulated coefficients of any real substance.
type stubModel struct {
	tsat float64
	err  error
}

func (s stubModel) LnPSat(float64) (float64, error)   { return 0, nil }
func (s stubModel) Pressure(float64) (float64, error) { return 0, nil }
func (s stubModel) ValidateTempRange(float64) bool    { return true }

func (s stubModel) Temperature(float64) (float64, error) {
	if s.err != nil {
		return 0, s.err
	}
	return s.tsat, nil
}

// TestInitialTemperatureGuessesBounds checks that the two guesses are the
// lowest and highest pure-component saturation temperatures.
//
// Under Raoult's law both the bubble and dew temperatures of a mixture
// lie between those of its pure components, so this pair brackets the
// root the secant solver is started on.
func TestInitialTemperatureGuessesBounds(t *testing.T) {
	testCases := []struct {
		name      string
		tsat      []float64
		low, high float64
	}{
		{"ascending", []float64{60, 80, 100}, 60, 100},
		{"descending", []float64{100, 80, 60}, 60, 100},
		{"unordered", []float64{80, 60, 100}, 60, 100},
		{"two components", []float64{56.2, 100}, 56.2, 100},
		{"negative temperatures", []float64{-40, -10}, -40, -10},
		{"extremes at the ends", []float64{100, 70, 80, 20}, 20, 100},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			models := make([]antoine.Model, len(tc.tsat))
			for i, tsat := range tc.tsat {
				models[i] = stubModel{tsat: tsat}
			}

			gotLow, gotHigh, err := InitialTemperatureGuesses(101.325, models)
			if err != nil {
				t.Fatalf("InitialTemperatureGuesses returned an unexpected error: %v", err)
			}

			if math.Abs(gotLow-tc.low) > 1e-12 {
				t.Errorf("lower guess = %g; want %g", gotLow, tc.low)
			}

			if math.Abs(gotHigh-tc.high) > 1e-12 {
				t.Errorf("upper guess = %g; want %g", gotHigh, tc.high)
			}

			if gotLow >= gotHigh {
				t.Errorf("guesses %g and %g are not distinct and ordered", gotLow, gotHigh)
			}
		})
	}
}

// TestInitialTemperatureGuessesIdenticalComponents checks that a mixture
// whose components boil at the same temperature is rejected.
//
// The secant method needs two distinct starting points; identical ones
// give a zero denominator on the first step.
func TestInitialTemperatureGuessesIdenticalComponents(t *testing.T) {
	testCases := []struct {
		name   string
		models []antoine.Model
	}{
		{
			name: "two components boiling alike",
			models: []antoine.Model{
				stubModel{tsat: 80},
				stubModel{tsat: 80},
			},
		},
		{
			// A single component leaves nothing to bracket between.
			name:   "one component",
			models: []antoine.Model{stubModel{tsat: 80}},
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, _, err := InitialTemperatureGuesses(101.325, tc.models); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestInitialTemperatureGuessesPropagatesError checks that a failure from
// a correlation reaches the caller rather than being absorbed into a
// guess of zero.
func TestInitialTemperatureGuessesPropagatesError(t *testing.T) {
	sentinel := errors.New("correlation failed")

	models := []antoine.Model{
		stubModel{tsat: 60},
		stubModel{err: sentinel},
	}

	_, _, err := InitialTemperatureGuesses(101.325, models)

	if !errors.Is(err, sentinel) {
		t.Errorf("error = %v; want the sentinel from the correlation", err)
	}
}

// TestInitialTemperatureGuessesMalformedInput checks the guards on the
// supplied correlations.
//
// The number of components is taken from the slice itself, so a stated
// count can no longer disagree with what is passed. That mismatch used
// to be representable and had two failure modes: indexing past the end
// of the working slice, or leaving it padded with zeros which were then
// reported as the lowest saturation temperature — a plausible-looking
// guess with nothing behind it.
func TestInitialTemperatureGuessesMalformedInput(t *testing.T) {
	testCases := []struct {
		name   string
		models []antoine.Model
	}{
		{"nil slice", nil},
		{"empty slice", []antoine.Model{}},
		{"nil model among valid ones", []antoine.Model{stubModel{tsat: 60}, nil}},
		{"only a nil model", []antoine.Model{nil}},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			if _, _, err := InitialTemperatureGuesses(101.325, tc.models); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestInitialTemperatureGuessesRealSubstances checks the guesses against
// two tabulated correlations, confirming the bracket is the one the VLE
// solvers expect.
//
// At one atmosphere a benzene/toluene mixture must boil somewhere
// between the two pure boiling points, so those are the bounds returned.
func TestInitialTemperatureGuessesRealSubstances(t *testing.T) {
	models := []antoine.Model{antoine.Benzene, antoine.Toluene}

	low, high, err := InitialTemperatureGuesses(101.325, models)
	if err != nil {
		t.Fatalf("InitialTemperatureGuesses returned an unexpected error: %v", err)
	}

	// Benzene is the more volatile component, so it sets the lower bound.
	if math.Abs(low-antoine.Benzene.Tn) > 0.5 {
		t.Errorf("lower guess = %.3f °C; want the benzene boiling point %.3f", low, antoine.Benzene.Tn)
	}

	if math.Abs(high-antoine.Toluene.Tn) > 0.5 {
		t.Errorf("upper guess = %.3f °C; want the toluene boiling point %.3f", high, antoine.Toluene.Tn)
	}

	if low >= high {
		t.Errorf("guesses %.3f and %.3f are not ordered", low, high)
	}
}
