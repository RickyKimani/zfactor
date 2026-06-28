package raoult

// PressureInput provides a composition vector and saturation pressures
// for Raoult-law pressure calculations.
//
// The composition represents:
//
//   - Liquid mole fractions x for bubble-pressure calculations.
//   - Vapor mole fractions y for dew-pressure calculations.
type PressureInput interface {
	Composition() []float64
	PSat() ([]float64, error)
}

// SaturationPressureInput uses user-supplied saturation pressures.
type SaturationPressureInput struct {
	// Compositions represent a composition vector. x for BUBL P calculations and y for DEW P calculations
	Compositions []float64
	PSats        []float64
}

// Composition returns the liquid composition.
func (p SaturationPressureInput) Composition() []float64 {
	return p.Compositions
}

// PSat returns the supplied saturation pressures.
func (p SaturationPressureInput) PSat() ([]float64, error) {
	return p.PSats, nil
}

// BubblePResult contains the bubble pressure and vapor composition.
type BubblePResult struct {
	P float64   // bubble pressure
	Y []float64 // vapor composition
}

// BubbleP calculates the bubble pressure and equilibrium vapor
// composition using Raoult's law.
//
//	P = Σ xi Pi_sat
//
//	yi = xi Pi_sat / P
func BubbleP(input PressureInput) (BubblePResult, error) {
	res, err := preparePressureInput(input)
	if err != nil {
		return BubblePResult{}, err
	}

	x := res.comp
	psat := res.psat
	n := res.n

	var p float64

	for i := range n {
		p += x[i] * psat[i]
	}

	y := make([]float64, n)

	for i := range n {
		y[i] = x[i] * psat[i] / p
	}

	return BubblePResult{
		P: p,
		Y: y,
	}, nil
}

// DewPResult contains the dew pressure and equilibrium liquid composition.
type DewPResult struct {
	P float64   // dew pressure
	X []float64 // liquid composition
}

// DewP calculates the dew pressure and equilibrium liquid composition
// using Raoult's law.
//
// The dew pressure is obtained from:
//
//	1/P = Σ (yi / Pi_sat)
//
// or equivalently:
//
//	P = 1 / Σ (yi / Pi_sat)
//
// where yi is the vapor-phase mole fraction and Pi_sat is the
// saturation pressure of component i at the specified temperature.
//
// The equilibrium liquid composition is then calculated from:
//
//	xi = yi P / Pi_sat
//
// The returned DewPResult contains the dew pressure and the liquid
// composition in equilibrium with the specified vapor composition.
func DewP(input PressureInput) (DewPResult, error) {
	res, err := preparePressureInput(input)
	if err != nil {
		return DewPResult{}, err
	}

	y := res.comp
	psat := res.psat
	n := res.n

	var denom float64
	for i := range n {
		denom += y[i] / psat[i]
	}

	p := 1 / denom

	x := make([]float64, n)
	for i := range n {
		x[i] = y[i] * p / psat[i]
	}

	return DewPResult{
		P: p, X: x,
	}, nil

}
