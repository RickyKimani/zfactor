// Package state provides functionality for defining thermodynamic states and generating
// visual representations such as PV diagrams.
package state

import (
	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/substance"
)

// State represents a specific thermodynamic state of a substance defined by its
// temperature and pressure.
type State struct {
	Substance   *substance.Substance
	Temperature float64 // Temperature in Kelvin
	Pressure    float64 // Pressure in bar
}

// NewState creates a new State object. It validates that the temperature and pressure
// are positive values.
func NewState(substance *substance.Substance, T, P float64) (*State, error) {
	if T <= 0 {
		return nil, zfactor.ErrTemp
	}
	if P <= 0 {
		return nil, zfactor.ErrPressure
	}
	return &State{
		Substance:   substance,
		Temperature: T,
		Pressure:    P,
	}, nil
}
