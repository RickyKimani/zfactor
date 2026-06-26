package raoult

import (
	"errors"

	"github.com/rickykimani/zfactor/antoine"
)

// TemperatureInput supplies the information required for bubble- and
// dew-temperature calculations using Raoult's law.
//
// The composition represents:
//
//   - Liquid mole fractions (x) for BUBL T.
//   - Vapor mole fractions (y) for DEW T.
type TemperatureInput interface {
	Composition() []float64
	Pressure() float64
	AntoineModels() []antoine.Model
	SolverOptions() SolverOptions
}

// saturationPressure returns the saturation saturationPressure of a component at temperature T.
//
// Temperatures outside the recommended Antoine correlation range return the
// computed saturation saturationPressure together with a *antoine.RangeError*. Since the
// correlation remains mathematically defined outside its fitted range, the
// saturationPressure is accepted and only non-range errors are propagated.
func saturationPressure(model antoine.Model, T float64) (float64, error) {
	psat, err := model.Pressure(T)

	var rerr *antoine.RangeError
	if err != nil && !errors.As(err, &rerr) {
		return 0, err
	}

	return psat, nil
}

// initialTemperatureGuesses returns two initial temperature estimates for the
// secant solver.
//
// The guesses are chosen as the minimum and maximum pure-component saturation
// temperatures corresponding to the specified system pressure. For Raoult's
// law, both the bubble and dew temperatures lie within these bounds.
func initialTemperatureGuesses(
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

// BubbleTResult contains the bubble temperature and the corresponding
// equilibrium vapor composition.
type BubbleTResult struct {
	T float64   // bubble temperature (°C)
	Y []float64 // equilibrium vapor composition
}

// bubbleResidual evaluates the Raoult's law bubble-temperature residual:
//
//	Σ xi Pi_sat(T) - P
//
// The bubble temperature is obtained by finding the temperature at which this
// residual is zero.
func bubbleResidual(
	T float64,
	P float64,
	x []float64,
	models []antoine.Model,
) (float64, error) {
	sum := 0.0

	for i, model := range models {
		psat, err := saturationPressure(model, T)
		if err != nil {
			return 0, err
		}
		sum += x[i] * psat
	}

	return sum - P, nil
}

// BubbleT calculates the bubble temperature and equilibrium vapor composition
// of a multicomponent mixture using Raoult's law.
//
// The bubble temperature is the temperature at which the first bubble of vapor
// forms from a liquid mixture at the specified pressure. It is obtained by
// solving
//
//	Σ xi Pi_sat(T) = P
//
// using the secant method.
func BubbleT(input TemperatureInput) (BubbleTResult, error) {
	x, p, models, opts, n, err := prepareTemperatureInput(input)
	if err != nil {
		return BubbleTResult{}, err
	}

	t0, t1, err := initialTemperatureGuesses(p, n, models)
	if err != nil {
		return BubbleTResult{}, err
	}

	T, err := secant(
		func(T float64) (float64, error) {
			return bubbleResidual(T, p, x, models)
		},
		t0,
		t1,
		opts,
	)

	if err != nil {
		return BubbleTResult{}, err
	}

	psat := make([]float64, n)
	for i, model := range models {
		psat[i], err = saturationPressure(model, T)
		if err != nil {
			return BubbleTResult{}, err
		}
	}

	y := make([]float64, n)
	for i := range n {
		y[i] = x[i] * psat[i] / p
	}

	return BubbleTResult{
		T: T,
		Y: y,
	}, nil
}

// DewTResult contains the dew temperature and the corresponding equilibrium
// liquid composition.
type DewTResult struct {
	T float64   // dew temperature (°C)
	X []float64 // equilibrium liquid composition
}

// dewResidual evaluates the Raoult's law dew-temperature residual:
//
//	Σ yi P / Pi_sat(T) - 1
//
// The dew temperature is obtained by finding the temperature at which this
// residual is zero.
func dewResidual(
	T float64,
	P float64,
	y []float64,
	models []antoine.Model,
) (float64, error) {
	sum := 0.0

	for i, model := range models {
		psat, err := saturationPressure(model, T)
		if err != nil {
			return 0, err
		}
		sum += y[i] * P / psat
	}
	return sum - 1, nil
}

// DewT calculates the dew temperature and equilibrium liquid composition
// of a multicomponent mixture using Raoult's law.
//
// The dew temperature is the temperature at which the first drop of liquid
// condenses from a vapor mixture at the specified pressure. It is obtained by
// solving
//
//	Σ yi P / Pi_sat(T) = 1
//
// using the secant method.
func DewT(input TemperatureInput) (DewTResult, error) {
	y, p, models, opts, n, err := prepareTemperatureInput(input)
	if err != nil {
		return DewTResult{}, err
	}

	t0, t1, err := initialTemperatureGuesses(p, n, models)
	if err != nil {
		return DewTResult{}, err
	}

	T, err := secant(
		func(T float64) (float64, error) {
			return dewResidual(T, p, y, models)
		},
		t0,
		t1,
		opts,
	)

	if err != nil {
		return DewTResult{}, err
	}

	psat := make([]float64, n)
	for i, model := range models {
		psat[i], err = saturationPressure(model, T)
		if err != nil {
			return DewTResult{}, err
		}
	}

	x := make([]float64, n)
	for i := range n {
		x[i] = y[i] * p / psat[i]
	}

	sum := 0.0
	for _, xi := range x {
		sum += xi
	}

	for i := range n {
		x[i] /= sum
	}

	return DewTResult{
		T: T,
		X: x,
	}, nil

}
