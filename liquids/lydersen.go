// Package liquids provides correlations and data for calculating liquid properties.
//
// Note: The Lydersen chart implementation relies on digitized data. Users should use
// the ReducedDensity function with care, as the underlying data values may change
// in subsequent versions as digitization accuracy improves or as better data sources
// are integrated.
package liquids

import (
	"fmt"
	"sort"
)

type point struct {
	Pr   float64
	RhoR float64
}

type isotherm struct {
	Tr     float64
	Points []point
}

type LydersenTable struct {
	Saturation []point
	Isotherms  []isotherm
}

// ReducedDensity calculates the reduced density (rho_r) for a given reduced temperature (Tr)
// and reduced pressure (Pr) using the Lydersen chart data.
// It performs bilinear interpolation between isotherms and pressure points.
func ReducedDensity(Tr, Pr float64) (float64, error) {
	isotherms := lydersenData.Isotherms
	if len(isotherms) == 0 {
		return 0, fmt.Errorf("lydersen table is empty")
	}

	// 1. Find the relevant isotherms (Tr interpolation)
	// Search for the first isotherm with Tr >= requested Tr
	idx := sort.Search(len(isotherms), func(i int) bool {
		return isotherms[i].Tr >= Tr
	})

	// Case: Tr is above the highest isotherm
	if idx == len(isotherms) {
		return 0, fmt.Errorf("Tr %g is above the maximum defined Tr (%g) in Lydersen table", Tr, isotherms[len(isotherms)-1].Tr)
	}

	// Case: Exact Tr match
	if isotherms[idx].Tr == Tr {
		return interpolatePr(isotherms[idx].Points, Pr)
	}

	// Case: Tr is below the lowest isotherm
	if idx == 0 {
		return 0, fmt.Errorf("Tr %g is below the minimum defined Tr (%g) in Lydersen table", Tr, isotherms[0].Tr)
	}

	// Case: Interpolate between two isotherms (idx-1 and idx)
	isoLow := isotherms[idx-1]
	isoHigh := isotherms[idx]

	rhoLow, err := interpolatePr(isoLow.Points, Pr)
	if err != nil {
		return 0, fmt.Errorf("failed to interpolate at lower Tr %g: %w", isoLow.Tr, err)
	}

	rhoHigh, err := interpolatePr(isoHigh.Points, Pr)
	if err != nil {
		return 0, fmt.Errorf("failed to interpolate at higher Tr %g: %w", isoHigh.Tr, err)
	}

	// Linear interpolation for Tr
	frac := (Tr - isoLow.Tr) / (isoHigh.Tr - isoLow.Tr)
	return rhoLow + frac*(rhoHigh-rhoLow), nil
}

// interpolatePr finds the density at a specific Pr within a single isotherm points slice
func interpolatePr(points []point, Pr float64) (float64, error) {
	if len(points) == 0 {
		return 0, fmt.Errorf("empty isotherm points")
	}

	// Search for the first point with Pr >= requested Pr
	idx := sort.Search(len(points), func(i int) bool {
		return points[i].Pr >= Pr
	})

	// Case: Pr is above the highest point
	if idx == len(points) {
		return 0, fmt.Errorf("Pr %g is above the maximum defined Pr (%g) for this isotherm", Pr, points[len(points)-1].Pr)
	}

	// Case: Exact Pr match
	if points[idx].Pr == Pr {
		return points[idx].RhoR, nil
	}

	// Case: Pr is below the lowest point
	if idx == 0 {
		return 0, fmt.Errorf("Pr %g is below the minimum defined Pr (%g) for this isotherm", Pr, points[0].Pr)
	}

	// Case: Interpolate between two points (idx-1 and idx)
	pLow := points[idx-1]
	pHigh := points[idx]

	// Linear interpolation for Pr
	frac := (Pr - pLow.Pr) / (pHigh.Pr - pLow.Pr)
	return pLow.RhoR + frac*(pHigh.RhoR-pLow.RhoR), nil
}
