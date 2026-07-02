package internal

import (
	"errors"

	"github.com/rickykimani/zfactor/antoine"
)

// InitialTemperatureGuesses returns two initial temperature estimates for the
// secant solver.
//
// The guesses are chosen as the minimum and maximum pure-component saturation
// temperatures corresponding to the specified system pressure. For Raoult's
// law, both the bubble and dew temperatures lie within these bounds.
func InitialTemperatureGuesses(
	P float64,
	n int,
	models []antoine.Model,
) (float64, float64, error) {
	var err error
	tsat := make([]float64, n)
	for i, model := range models {
		tsat[i], err = model.Temperature(P)
		if err != nil {
			return 0, 0, err
		}
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
