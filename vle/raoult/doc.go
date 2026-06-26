// Package raoult implements vapor-liquid equilibrium (VLE) calculations
// based on Raoult's law for ideal liquid mixtures.
//
// The package provides routines for:
//
//   - Bubble-pressure (BubbleP)
//   - Dew-pressure (DewP)
//   - Bubble-temperature (BubbleT)
//   - Dew-temperature (DewT)
//
// Pressure-based calculations may use either user-supplied saturation
// pressures or Antoine vapor-pressure correlations. Temperature-based
// calculations use Antoine correlations together with a numerical root
// solver to determine equilibrium temperatures.
//
// These routines assume:
//
//   - An ideal liquid phase.
//
//   - An ideal vapor phase.
//
//   - Phase equilibrium described by Raoult's law:
//
//     yᵢ P = xᵢ Pᵢˢᵃᵗ
//
// For non-ideal liquid mixtures, see the activity coefficient models in
// the activity package and the modified Raoult's law implementations.
package raoult
