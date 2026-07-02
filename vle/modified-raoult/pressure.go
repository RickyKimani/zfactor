package modified_raoult

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor/activity"
	"github.com/rickykimani/zfactor/vle"
)

// PressureInput supplies the information required for bubble- and
// dew-pressure calculations using the modified Raoult's law.
//
// The composition represents:
//
//   - Liquid mole fractions (x) for BUBL P.
//   - Vapor mole fractions (y) for DEW P.
type PressureInput interface {
	Temperature() float64
	Composition() []float64
	PSat() ([]float64, error)
	ActivityModel() activity.Model
	SolverOptions() vle.SolverOptions
}

// SaturationPressureInput uses user-supplied saturation pressures.
type SaturationPressureInput struct {

	// Temperature (°C).
	T float64

	// Composition vector.
	Compositions []float64

	// Saturation pressures (kPa).
	PSats []float64

	// Activity-coefficient model.
	Activity activity.Model

	Options vle.SolverOptions
}

func (s SaturationPressureInput) Temperature() float64 {
	return s.T
}

func (s SaturationPressureInput) Composition() []float64 {
	return s.Compositions
}

func (s SaturationPressureInput) PSat() ([]float64, error) {
	return s.PSats, nil
}

func (s SaturationPressureInput) ActivityModel() activity.Model {
	return s.Activity
}

func (s SaturationPressureInput) SolverOptions() vle.SolverOptions {
	return s.Options
}

// BubblePResult contains the bubble pressure and equilibrium vapor
// composition.
type BubblePResult struct {
	P float64   // bubble pressure (kPa)
	Y []float64 // equilibrium vapor composition
}

// BubbleP calculates the bubble pressure and equilibrium vapor
// composition using the modified Raoult's law.
//
// The bubble pressure is obtained from
//
//	P = Σ xi γi Pi_sat
//
// where γi is evaluated from the supplied activity model at the
// specified liquid composition and temperature.
//
// The equilibrium vapor composition is then
//
//	yi = xi γi Pi_sat / P.
func BubbleP(input PressureInput) (BubblePResult, error) {
	res, err := preparePressureInput(input)
	if err != nil {
		return BubblePResult{}, err
	}

	gamma, err := res.model.Activity()
	if err != nil {
		return BubblePResult{}, err
	}

	var P float64
	for i := range res.n {
		P += res.comp[i] * gamma[i] * res.psat[i]
	}

	y := make([]float64, res.n)

	for i := range res.n {
		y[i] = res.comp[i] * gamma[i] * res.psat[i] / P
	}

	return BubblePResult{
		P: P,
		Y: y,
	}, nil
}

// DewPResult contains the dew pressure and equilibrium liquid
// composition.
type DewPResult struct {
	P float64   // dew pressure (kPa)
	X []float64 // equilibrium liquid composition
}

// DewP calculates the dew pressure and equilibrium liquid
// composition using the modified Raoult's law.
func DewP(input PressureInput) (DewPResult, error) {
	res, err := preparePressureInput(input)
	if err != nil {
		return DewPResult{}, err
	}

	y := res.comp
	psat := res.psat
	model := res.model
	n := res.n

	gamma := make([]float64, n)
	for i := range gamma {
		gamma[i] = 1
	}

	x := make([]float64, n)

	tol := res.opts.Tol()
	maxIter := res.opts.MaxIter()

	for iter := 0; iter < maxIter; iter++ {

		// Compute pressure.
		var denom float64
		for i := range n {
			denom += y[i] / (gamma[i] * psat[i])
		}

		P := 1 / denom

		// Compute liquid composition.
		var sum float64
		for i := range n {
			x[i] = y[i] * P / (gamma[i] * psat[i])
			sum += x[i]
		}

		// Normalize.
		for i := range n {
			x[i] /= sum
		}

		// Evaluate activity coefficients.
		model = model.WithComposition(x)

		gammaNew, err := model.Activity()
		if err != nil {
			return DewPResult{}, err
		}

		var maxDiff float64
		for i := range n {
			diff := math.Abs(gammaNew[i] - gamma[i])
			if diff > maxDiff {
				maxDiff = diff
			}
		}

		gamma = gammaNew

		if maxDiff < tol {
			return DewPResult{
				P: P,
				X: x,
			}, nil
		}
	}

	return DewPResult{}, errors.New(
		"modified Raoult dew-pressure calculation failed to converge",
	)
}
