// Package cp provides coefficients for the calculation of heat capacities
// of gases (in the ideal-gas state), solids, and liquids.
//
// The parameters stored in HeatCapacity struct correspond to the equation:
//
//	Cp/R = A + B*T + C*T^2 + D*T^-2
//
// Note: This package currently only provides the data constants. Calculation
// functions for integrals (Enthalpy/Entropy changes) are pending implementation.
package cp

// HeatCapacity holds the constants for the heat capacity equation.
// Cp/R = A + B*T + C*T^2 + D*T^-2
type HeatCapacity struct {
	Name    string
	Formula string
	TMin    float64
	TMax    float64
	Cp298   float64 // Cp/R at 298.15 K
	A       float64
	B       float64
	C       float64
	D       float64
}
