package modified_raoult

import (
	"errors"

	"github.com/rickykimani/zfactor/activity"
	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/vle"
	"github.com/rickykimani/zfactor/vle/internal"
)

// TemperatureInput supplies the information required for bubble- and
// dew-temperature calculations using modified Raoult's law.
//
// The composition represents:
//
//   - Liquid mole fractions (x) for BUBL T.
//   - Vapor mole fractions (y) for DEW T.
type TemperatureInput interface {
	Composition() []float64
	Pressure() float64
	AntoineModels() []antoine.Model
	ActivityModel() activity.Model
	SolverOptions() vle.SolverOptions
}

// saturationPressure returns the saturation pressure of a component at temperature T.
//
// Temperatures outside the recommended Antoine correlation range return the
// computed saturation saturationPressure together with a *antoine.RangeError*. Since the
// correlation remains mathematically defined outside its fitted range, the
// saturation pressure is accepted and only non-range errors are propagated.
func saturationPressure(model antoine.Model, T float64) (float64, error) {
	psat, err := model.Pressure(T)

	var rerr *antoine.RangeError
	if err != nil && !errors.As(err, &rerr) {
		return 0, err
	}

	return psat, nil
}

// BubbleTResult contains the bubble temperature and the corresponding
// equilibrium vapor composition.
type BubbleTResult struct {
	T float64   // bubble temperature (°C)
	Y []float64 // equilibrium vapor composition
}

// bubbleResidual evaluates the modified Raoult's law bubble-temperature
// residual:
//
//	Σ xi γi(T,x) Pi_sat(T) - P
func bubbleResidual(
	T float64,
	P float64,
	models []antoine.Model,
	model activity.Model,
) (float64, error) {

	model = model.WithTemperature(toKelvin(T))

	gamma, err := model.Activity()
	if err != nil {
		return 0, err
	}

	x := model.Composition()

	var sum float64
	for i, am := range models {
		psat, err := saturationPressure(am, T)
		if err != nil {
			return 0, err
		}

		sum += x[i] * gamma[i] * psat
	}

	return sum - P, nil
}

// BubbleT calculates the bubble temperature and equilibrium vapor
// composition using the modified Raoult's law.
//
// The bubble temperature is obtained by solving
//
//	Σ xi γi(T,x) Pi_sat(T) = P
//
// using the secant method.
func BubbleT(input TemperatureInput) (BubbleTResult, error) {

	res, err := prepareTemperatureInput(input)
	if err != nil {
		return BubbleTResult{}, err
	}

	x := res.comp
	p := res.p
	n := res.n
	models := res.models // antoine
	activityModel := res.activity.WithComposition(x)
	opts := res.opts

	t0, t1, err := internal.InitialTemperatureGuesses(p, n, models)
	if err != nil {
		return BubbleTResult{}, err
	}

	T, err := internal.Secant(
		func(T float64) (float64, error) {
			return bubbleResidual(
				T,
				p,
				models,
				activityModel,
			)
		},
		t0,
		t1,
		opts,
	)
	if err != nil {
		return BubbleTResult{}, err
	}

	// Evaluate γ at the converged temperature.
	activityModel = activityModel.WithTemperature(toKelvin(T))

	gamma, err := activityModel.Activity()
	if err != nil {
		return BubbleTResult{}, err
	}

	psat := make([]float64, n)
	for i, am := range models {
		psat[i], err = saturationPressure(am, T)
		if err != nil {
			return BubbleTResult{}, err
		}
	}

	y := make([]float64, n)
	for i := range n {
		y[i] = x[i] * gamma[i] * psat[i] / p
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

// dewLiquidComposition performs the inner fixed-point iteration required
// for modified Raoult's law dew-temperature calculations.
//
// For a specified temperature and pressure, the liquid composition is
// obtained by iterating between
//
//	xi = yi P / (γi Pi_sat)
//
// and evaluation of the activity coefficients until convergence.
//
// The returned residual is
//
//	Σ yi P / (γi Pi_sat) - 1
//
// which is used by the outer secant solver.
func dewLiquidComposition(
	T float64,
	P float64,
	y []float64,
	models []antoine.Model,
	activityModel activity.Model,
	opts vle.SolverOptions,
) (
	x []float64,
	gamma []float64,
	residual float64,
	err error,
) {
	n := len(y)

	psat := make([]float64, n)
	for i, am := range models {
		psat[i], err = saturationPressure(am, T)
		if err != nil {
			return
		}
	}

	// Initial guess assuming γ = 1.
	x0 := make([]float64, n)
	sum := 0.0

	for i := range n {
		x0[i] = y[i] * P / psat[i]
		sum += x0[i]
	}

	for i := range n {
		x0[i] /= sum
	}

	x, err = internal.FixedPoint(
		x0,
		func(x []float64) ([]float64, error) {

			m := activityModel.
				WithTemperature(toKelvin(T)).
				WithComposition(x)

			gamma, err = m.Activity()
			if err != nil {
				return nil, err
			}

			next := make([]float64, n)

			sum := 0.0
			for i := range n {
				next[i] = y[i] * P /
					(gamma[i] * psat[i])
				sum += next[i]
			}

			for i := range n {
				next[i] /= sum
			}

			return next, nil
		},
		opts,
	)
	if err != nil {
		return
	}

	activityModel = activityModel.
		WithTemperature(toKelvin(T)).
		WithComposition(x)

	gamma, err = activityModel.Activity()
	if err != nil {
		return
	}

	residual = 0.0
	for i := range n {
		residual += y[i] * P /
			(gamma[i] * psat[i])
	}
	residual--

	return
}

// dewResidual evaluates the modified Raoult's law dew-temperature residual:
//
//	Σ yi P / (γi(T,x) Pi_sat(T)) - 1
func dewResidual(
	T float64,
	P float64,
	y []float64,
	models []antoine.Model,
	activityModel activity.Model,
	opts vle.SolverOptions,
) (float64, error) {

	_, _, residual, err := dewLiquidComposition(
		T,
		P,
		y,
		models,
		activityModel,
		opts,
	)

	return residual, err
}

// DewT calculates the dew temperature and equilibrium liquid composition
// using the modified Raoult's law.
//
// The dew temperature is obtained by solving
//
//	Σ yi P / (γi(T,x) Pi_sat(T)) = 1
//
// using the secant method. Since the activity coefficients depend on the
// unknown liquid composition, each residual evaluation performs an inner
// fixed-point iteration.
func DewT(input TemperatureInput) (DewTResult, error) {

	res, err := prepareTemperatureInput(input)
	if err != nil {
		return DewTResult{}, err
	}

	y := res.comp
	p := res.p
	n := res.n
	models := res.models
	activityModel := res.activity
	opts := res.opts

	t0, t1, err := internal.InitialTemperatureGuesses(p, n, models)
	if err != nil {
		return DewTResult{}, err
	}

	T, err := internal.Secant(
		func(T float64) (float64, error) {
			return dewResidual(
				T,
				p,
				y,
				models,
				activityModel,
				opts,
			)
		},
		t0,
		t1,
		opts,
	)
	if err != nil {
		return DewTResult{}, err
	}

	x, _, _, err := dewLiquidComposition(
		T,
		p,
		y,
		models,
		activityModel,
		opts,
	)
	if err != nil {
		return DewTResult{}, err
	}

	return DewTResult{
		T: T,
		X: x,
	}, nil
}
