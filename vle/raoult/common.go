package raoult

import (
	"errors"
	"fmt"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/antoine"
)

const tolerance = 1e-6

// MixtureInput contains thermodynamic information available to a solver for performing saturation calculations
// - T is the temperature in °C (mostly/if the correlation originates from this library)
// and is ignored when performing Temperature calculations (BUBL T and DEW T)
//
// - P is the pressure in kPa (mostly/if the correlation originates from this library)
// and is ignored when performing Pressure calculations (BUBL P and DEW P)
type MixtureInput struct {
	T float64
	P float64

	// Compositions represent a composition vector. x for BUBL P calculations and y for DEW P calculations
	Compositions []float64

	Antoine []antoine.Model

	Options SolverOptions
}

// Composition returns the liquid composition.
func (m MixtureInput) Composition() []float64 {
	return m.Compositions
}

// PSat computes saturation pressures using Antoine correlations.
func (m MixtureInput) PSat() ([]float64, error) {
	if m.T <= -273.15 {
		return nil, zfactor.ErrTemp
	}

	n := len(m.Antoine)
	psat := make([]float64, n)

	for i, model := range m.Antoine {
		p, err := model.Pressure(m.T)
		if err != nil {
			return nil, err
		}
		psat[i] = p
	}

	return psat, nil
}

func (m MixtureInput) Pressure() float64 {
	return m.P
}

func (m MixtureInput) AntoineModels() []antoine.Model {
	return m.Antoine
}

func (m MixtureInput) SolverOptions() SolverOptions {
	return m.Options
}

// validateComposition validates a composition vector.
func validateComposition(w []float64) error {
	if len(w) == 0 {
		return errors.New("no components provided")
	}

	var sum float64

	for _, xi := range w {
		if xi < 0 || xi > 1 {
			return zfactor.ErrMolFracVal
		}
		sum += xi
	}

	if math.Abs(sum-1.0) > tolerance {
		return zfactor.ErrMolFracSum
	}

	return nil
}

// validatePSat validates saturation pressures.
func validatePSat(psat []float64) error {
	if len(psat) == 0 {
		return errors.New("no saturation pressures provided")
	}

	for _, p := range psat {
		if p <= 0 {
			return zfactor.ErrPressure
		}
	}

	return nil
}

func preparePressureInput(input PressureInput) (
	comp,
	psat []float64,
	n int,
	err error,
) {
	comp = input.Composition()
	if err = validateComposition(comp); err != nil {
		return nil, nil, 0, err
	}

	psat, err = input.PSat()
	if err != nil {
		return nil, nil, 0, err
	}

	n = len(comp)

	if n < 2 {
		return nil, nil, 0, errors.New(
			"pressure calculations require at least two components",
		)
	}

	if len(psat) != n {
		return nil, nil, 0, errors.New(
			"number of saturation pressures must match number of components",
		)
	}

	if err = validatePSat(psat); err != nil {
		return nil, nil, 0, err
	}

	return comp, psat, n, nil
}

func prepareTemperatureInput(
	input TemperatureInput,
) (
	comp []float64,
	P float64,
	models []antoine.Model,
	opts SolverOptions,
	n int,
	err error,
) {
	comp = input.Composition()

	if err = validateComposition(comp); err != nil {
		return nil, 0, nil, SolverOptions{}, 0, err
	}

	P = input.Pressure()
	if P <= 0 {
		return nil, 0, nil, SolverOptions{}, 0, zfactor.ErrPressure
	}

	models = input.AntoineModels()

	n = len(comp)

	if n < 2 {
		return nil, 0, nil, SolverOptions{}, 0, errors.New(
			"temperature calculations require at least two components",
		)
	}

	if len(models) != n {
		return nil, 0, nil, SolverOptions{}, 0, errors.New(
			"number of Antoine models must match number of components",
		)
	}

	for i, model := range models {
		if model == nil {
			return nil, 0, nil, SolverOptions{}, 0, fmt.Errorf(
				"antoine model %d is nil",
				i,
			)
		}
	}

	return comp, P, models, input.SolverOptions(), n, nil
}
