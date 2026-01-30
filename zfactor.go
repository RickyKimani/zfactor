// Copyright (c) 2025 Ricky Kimani
// Licensed under the MIT License. See LICENSE file in the project root for full license information.

// Package zfactor provides a comprehensive library for thermodynamic property calculations
// and visualization. It includes tools for solving Cubic Equations of State (EOS),
// estimating properties using correlations like Lee-Kesler, calculating liquid properties,
// and generating Pressure-Volume (PV) diagrams.
package zfactor

const (
	// RSI is the Universal Gas Constant in SI units [J/(mol·K)].
	RSI = 8.314

	// AtmPa is the standard atmospheric pressure in Pascals (Pa).
	AtmPa = 101_325.0

	// AtmKPa is the standard atmospheric pressure in Kilopascals (kPa).
	AtmKPa = AtmPa * 1e-3

	// AtmBar is the standard atmospheric pressure in Bars.
	AtmBar = AtmPa * 1e-5
)

// Args holds the thermodynamic state arguments to prevent order-dependent errors.
// It is used to pass parameters like Temperature and Pressure safely.
type Args struct {
	T float64 // Temperature
	P float64 // Pressure
	R float64 // Gas constant
	B float64 // Second virial coefficient
	C float64 // Third virial coefficient
}
