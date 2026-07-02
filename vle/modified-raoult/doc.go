// Package modified_raoult implements vapor-liquid equilibrium (VLE)
// calculations using the modified Raoult's law.
//
// Modified Raoult's law extends the ideal-solution assumption by
// accounting for liquid-phase non-ideality through activity
// coefficients:
//
//	yi P = xi γi Pi,sat
//
// where yi and xi are the vapor- and liquid-phase mole fractions,
// γi is the activity coefficient of component i, P is the system
// pressure, and Pi,sat is the saturation pressure of the pure
// component.
//
// Saturation pressures are obtained from Antoine correlations, while
// activity coefficients are supplied by an activity-coefficient model
// such as Wilson or NRTL.
//
// Temperatures supplied to this package are expressed in degrees
// Celsius to remain consistent with the Antoine package. Activity
// models internally receive temperatures in kelvin.
//
// The package provides routines for bubble- and dew-point pressure
// and temperature calculations for non-ideal liquid mixtures.
package modified_raoult
