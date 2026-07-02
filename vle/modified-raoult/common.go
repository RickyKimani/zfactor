package modified_raoult

import (
	"errors"
	"fmt"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/activity"
	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/vle"
)

const (
	tempConv  = 273.15
	tolerance = 1e-6
)

func toKelvin(tC float64) float64 {
	return tC + tempConv
}

// MixtureInput contains the thermodynamic information required for
// modified-Raoult VLE calculations.
//
// Temperatures are expressed in °C and pressures in kPa. Activity
// coefficient models are internally evaluated in Kelvin.
type MixtureInput struct {
	T float64
	P float64

	// x for bubble calculations, y for dew calculations.
	Compositions []float64

	// Pure-component vapor pressure correlations.
	Antoine []antoine.Model

	// Liquid-phase activity coefficient model.
	Activity activity.Model

	Options vle.SolverOptions
}

func (m MixtureInput) Temperature() float64 {
	return m.T
}

func (m MixtureInput) Pressure() float64 {
	return m.P
}

func (m MixtureInput) Composition() []float64 {
	return m.Compositions
}

func (m MixtureInput) ActivityModel() activity.Model {
	return m.Activity
}

func (m MixtureInput) AntoineModels() []antoine.Model {
	return m.Antoine
}

func (m MixtureInput) SolverOptions() vle.SolverOptions {
	return m.Options
}

// PSat evaluates the saturation pressures using the stored Antoine
// correlations.
func (m MixtureInput) PSat() ([]float64, error) {
	if m.T <= -tempConv {
		return nil, zfactor.ErrTemp
	}

	psat := make([]float64, len(m.Antoine))

	for i, model := range m.Antoine {
		p, err := saturationPressure(model, m.T)
		if err != nil {
			return nil, err
		}
		psat[i] = p
	}

	return psat, nil
}

// activityCoefficients evaluates the activity coefficients for a
// specified liquid composition and temperature.
//
// The supplied temperature is in °C and is converted internally to K.
func activityCoefficients(
	model activity.Model,
	T float64,
	x []float64,
) ([]float64, error) {

	if model == nil {
		return nil, errors.New("activity model is nil")
	}

	model = model.
		WithTemperature(toKelvin(T)).
		WithComposition(x)

	gamma, err := model.Activity()
	if err != nil {
		return nil, err
	}

	if err := validateGamma(gamma); err != nil {
		return nil, err
	}

	return gamma, nil
}

// validateComposition validates a composition vector.
func validateComposition(x []float64) error {

	if len(x) == 0 {
		return errors.New("no components provided")
	}

	sum := 0.0

	for _, xi := range x {
		if xi < 0 || xi > 1 {
			return zfactor.ErrMolFracVal
		}
		sum += xi
	}

	if math.Abs(sum-1) > tolerance {
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

// validateGamma validates activity coefficients.
func validateGamma(gamma []float64) error {

	if len(gamma) == 0 {
		return errors.New("no activity coefficients provided")
	}

	for _, g := range gamma {
		if g <= 0 {
			return errors.New("activity coefficients must be positive")
		}
	}

	return nil
}

type presPrepResult struct {
	T     float64
	comp  []float64
	psat  []float64
	model activity.Model
	n     int
	opts  vle.SolverOptions
}

func preparePressureInput(input PressureInput) (presPrepResult, error) {

	comp := input.Composition()
	if err := validateComposition(comp); err != nil {
		return presPrepResult{}, err
	}

	T := input.Temperature()
	if T <= -tempConv {
		return presPrepResult{}, zfactor.ErrTemp
	}

	psat, err := input.PSat()
	if err != nil {
		return presPrepResult{}, err
	}

	n := len(comp)

	if n < 2 {
		return presPrepResult{}, errors.New(
			"pressure calculations require at least two components",
		)
	}

	if len(psat) != n {
		return presPrepResult{}, errors.New(
			"number of saturation pressures must match number of components",
		)
	}

	if err := validatePSat(psat); err != nil {
		return presPrepResult{}, err
	}

	model := input.ActivityModel()
	if model == nil {
		return presPrepResult{}, errors.New("activity model is nil")
	}

	model = model.WithTemperature(toKelvin(T)).WithComposition(comp)
	return presPrepResult{
		T:     T,
		comp:  comp,
		psat:  psat,
		model: model,
		n:     n,
		opts:  input.SolverOptions(),
	}, nil
}

type tempPrepResult struct {
	comp     []float64
	p        float64
	models   []antoine.Model
	activity activity.Model
	opts     vle.SolverOptions
	n        int
}

func prepareTemperatureInput(input TemperatureInput) (tempPrepResult, error) {

	comp := input.Composition()
	if err := validateComposition(comp); err != nil {
		return tempPrepResult{}, err
	}

	p := input.Pressure()
	if p <= 0 {
		return tempPrepResult{}, zfactor.ErrPressure
	}

	models := input.AntoineModels()

	n := len(comp)

	if n < 2 {
		return tempPrepResult{}, errors.New(
			"temperature calculations require at least two components",
		)
	}

	if len(models) != n {
		return tempPrepResult{}, errors.New(
			"number of Antoine models must match number of components",
		)
	}

	for i, model := range models {
		if model == nil {
			return tempPrepResult{},
				fmt.Errorf("antoine model %d is nil", i)
		}
	}

	if input.ActivityModel() == nil {
		return tempPrepResult{}, errors.New(
			"activity model is nil",
		)
	}

	return tempPrepResult{
		comp:     comp,
		p:        p,
		models:   models,
		activity: input.ActivityModel(),
		opts:     input.SolverOptions(),
		n:        n,
	}, nil
}
