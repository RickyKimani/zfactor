package internal

import (
	"errors"
	"fmt"

	"github.com/rickykimani/zfactor/antoine"
)

// InitialTemperatureGuesses returns two initial temperature estimates for the
// secant solver.
//
// The guesses are chosen as the minimum and maximum pure-component saturation
// temperatures corresponding to the specified system pressure. For Raoult's
// law, both the bubble and dew temperatures lie within these bounds.
//
// The number of components is taken from models, so a caller cannot state a
// count that disagrees with the correlations it supplies.
//
// It returns an error if no models are given, if any is nil, if a correlation
// fails, or if every component shares a saturation temperature, since the
// secant method needs two distinct starting points.
func InitialTemperatureGuesses(
	P float64,
	models []antoine.Model,
) (float64, float64, error) {
	if len(models) == 0 {
		return 0, 0, errors.New("no components provided")
	}

	tsat := make([]float64, len(models))

	for i, model := range models {
		if model == nil {
			return 0, 0, fmt.Errorf("antoine model %d is nil", i)
		}

		t, err := model.Temperature(P)
		if err != nil {
			return 0, 0, err
		}

		tsat[i] = t
	}

	t0 := tsat[0]
	t1 := tsat[0]

	for _, t := range tsat[1:] {
		if t < t0 {
			t0 = t
		}
		if t > t1 {
			t1 = t
		}
	}

	if t0 == t1 {
		return 0, 0, errors.New(
			"unable to generate distinct initial guesses",
		)
	}

	return t0, t1, nil
}
