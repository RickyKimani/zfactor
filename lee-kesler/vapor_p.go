package leekesler

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor"
)

// normalErr is returned when the normal boiling point is invalid.
var normalErr error = errors.New("normal boiling point cannot be less than or equal to 0")

// lnReducedVaporPressureSimple calculates ln(Pr0) for the simple fluid term.
// Lee & Kesler (1975).
func lnReducedVaporPressureSimple(Tr float64) float64 {
	return 5.92714 - 6.09648/Tr - 1.28862*math.Log(Tr) + 0.169347*math.Pow(Tr, 6)
}

// lnReducedVaporPressureCorrection calculates ln(Pr1) for the correction fluid term.
// Lee & Kesler (1975).
func lnReducedVaporPressureCorrection(Tr float64) float64 {
	return 15.2518 - 15.6875/Tr - 13.4721*math.Log(Tr) + 0.43577*math.Pow(Tr, 6)
}

// EstimateAcentricFactor calculates the acentric factor (ω) using the Lee-Kesler correlation.
// It requires the Normal Boiling Point (Tn), Critical Temperature (Tc), and Critical Pressure (Pc).
//
// Arguments:
//   - Tn: Normal Boiling Point (K)
//   - Tc: Critical Temperature (K)
//   - Pc: Critical Pressure (bar) - Must be in bar as it relies on P_atm = 1.01325 bar.
//
// Formula: ω = (ln(Patm/Pc) - ln(Pr0_Tn)) / ln(Pr1_Tn)
func EstimateAcentricFactor(Tn, Tc, Pc float64) (float64, error) {
	if Tc <= 0 || Pc <= 0 {
		return 0, zfactor.ErrCriticalProp
	}
	if Tn <= 0 {
		return 0, normalErr
	}

	// Reduced normal boiling point
	Trn := Tn / Tc

	// ln(Pr_sat) at Normal Boiling Point. Pr = Patm / Pc.
	lnPrnSat := math.Log(zfactor.AtmBar / Pc)

	lnPr0 := lnReducedVaporPressureSimple(Trn)
	lnPr1 := lnReducedVaporPressureCorrection(Trn)

	// Avoid division by zero if lnPr1 is very small (unlikely for valid Trn < 1)
	if math.Abs(lnPr1) < 1e-9 {
		return 0, errors.New("calculation error: reference fluid vapor pressure term is zero")
	}

	return (lnPrnSat - lnPr0) / lnPr1, nil
}

// VaporPressure estimates the saturation vapor pressure (Psat) at a given temperature T.
// This function internally calculates the acentric factor using the Normal Boiling Point (Tn)
// to ensure consistency with the Lee-Kesler correlation.
//
// Arguments:
//   - T: Temperature (K)
//   - Tn: Normal Boiling Point (K)
//   - Tc: Critical Temperature (K)
//   - Pc: Critical Pressure (bar) - Must be in bar for correct acentric factor estimation.
//
// Returns Psat in bar.
func VaporPressure(T, Tn, Tc, Pc float64) (float64, error) {
	if T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if Tc <= 0 || Pc <= 0 {
		return 0, zfactor.ErrCriticalProp
	}

	Tr := T / Tc

	omega, err := EstimateAcentricFactor(Tn, Tc, Pc)
	if err != nil {
		return 0, err
	}

	lnPr0 := lnReducedVaporPressureSimple(Tr)
	lnPr1 := lnReducedVaporPressureCorrection(Tr)

	// ln(Pr) = ln(Pr0) + omega * ln(Pr1)
	lnPr := lnPr0 + omega*lnPr1

	// Psat = Pc * exp(lnPr)
	return Pc * math.Exp(lnPr), nil
}
