// Package ov provides the Orbey-Vera generalized correlation for the third
// virial coefficient.
//
// The reduced third virial coefficient is calculated as:
//
//	C * Pc^2 / (R^2 * Tc^2) = C0 + omega * C1
//
// where C0 is the simple-fluid contribution and C1 is the correction for the
// acentric factor (omega). The second virial coefficient in the same reduced
// formulation is supplied by the Abbott package.
//
// The residual-property functions use the Leiden form of the virial equation,
// truncated after the third virial coefficient. Given reduced temperature Tr
// and pressure Pr, Delta solves for the smallest positive reduced-density root
// from:
//
//	C * Delta^3 + B * Delta^2 + Delta - Pr/Tr = 0
//
// This is the vapor-like root described by the virial formulation. The
// residual enthalpy returned by ResidualEnthalpy is H^R/(R*Tc), and the
// residual entropy returned by ResidualEntropy is S^R/R.
package ov

import (
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/abbott"
)

// C0 calculates the simple-fluid contribution to the reduced third virial
// coefficient.
//
//	C0 = 0.01407 + 0.02432/Tr - 0.00313/Tr^10.5
//
// It returns an error if Tr <= 0.
func C0(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return 0.01407 + 0.02432/Tr - 0.00313/math.Pow(Tr, 10.5), nil
}

// C1 calculates the acentric-factor correction to the reduced third virial
// coefficient.
//
//	C1 = -0.02676 + 0.0177/Tr^2.8 + 0.04/Tr^3
//	     - 0.003/Tr^6 - 0.00228/Tr^10.5
//
// It returns an error if Tr <= 0.
func C1(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return -0.02676 + 0.0177/math.Pow(Tr, 2.8) + 0.04/math.Pow(Tr, 3) - 0.003/math.Pow(Tr, 6) - 0.00228/math.Pow(Tr, 10.5), nil
}

// ReducedC calculates the reduced third virial coefficient.
//
//	C * Pc^2 / (R^2 * Tc^2) = C0 + omega*C1
//
// It returns an error if Tr <= 0.
func ReducedC(Tr, acentric float64) (float64, error) {
	// No need for checking Tr validity here
	C0, err := C0(Tr)
	if err != nil {
		return 0, err
	}
	C1, err := C1(Tr)
	if err != nil {
		return 0, err
	}
	return C0 + acentric*C1, nil
}

// DC0 calculates the derivative of C0 with respect to reduced temperature.
//
//	dC0/dTr = -0.02432/Tr^2 + 0.032865/Tr^11.5
//
// It returns an error if Tr <= 0.
func DC0(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return -0.02432/math.Pow(Tr, 2) + 10.5*0.00313/math.Pow(Tr, 11.5), nil
}

// DC1 calculates the derivative of C1 with respect to reduced temperature.
//
//	dC1/dTr = -0.04956/Tr^3.8 - 0.12/Tr^4 + 0.018/Tr^7
//	           + 0.02394/Tr^11.5
//
// It returns an error if Tr <= 0.
func DC1(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return -2.8*0.0177/math.Pow(Tr, 3.8) - 3*0.04/math.Pow(Tr, 4) + 6*0.003/math.Pow(Tr, 7) + 10.5*0.00228/math.Pow(Tr, 11.5), nil
}

// DeltaArgs holds the reduced data necessary to calculate Delta
type DeltaArgs struct {
	// C is the reduced third virial coefficient. Use [ReducedC] to calculate it.
	C float64
	// B is the reduced second virial coefficient. Use [abbott.ReducedB] to calculate it.
	B float64
	// Pr is the reduced pressure Pr = P/Pc
	Pr float64
	// Tr is the reduced temperature Tr = T/Tc
	Tr float64
}

// Delta solves the truncated Leiden virial equation for the smallest positive
// reduced-density root:
//
//	C*Delta^3 + B*Delta^2 + Delta - Pr/Tr = 0
//
// Complex and non-positive roots are discarded. An error is returned if the
// reduced pressure or temperature is not positive, if the cubic solver fails,
// or if no physically meaningful root remains.
func Delta(d DeltaArgs) (float64, error) {
	if d.Pr <= 0 {
		return 0, zfactor.ErrInvalidPr
	}
	if d.Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	deltas, err := zfactor.SolveCubic(d.C, d.B, 1, -d.Pr/d.Tr)
	if err != nil {
		return 0, err
	}
	delta, err := clean(deltas)
	if err != nil {
		return 0, err
	}
	return delta, nil
}

// ResidualEnthalpy calculates the dimensionless residual enthalpy H^R/(R*Tc)
// using the second- and third-virial correlations.
//
// The result is evaluated as:
//
//	H^R/(R*Tc) = Tr * (Delta*B_H + Delta^2*C_H)
//
// It returns an error if Tr or Pr is not positive or if the reduced-density
// equation has no physically meaningful solution.
func ResidualEnthalpy(Tr, Pr, acentric float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}
	redC, err := ReducedC(Tr, acentric)
	if err != nil {
		return 0, nil
	}
	redB, err := abbott.ReducedB(Tr, acentric)
	if err != nil {
		return 0, err
	}
	BH, err := bH(Tr, acentric)
	if err != nil {
		return 0, err
	}
	CH, err := cH(Tr, acentric)

	delta, err := Delta(DeltaArgs{
		C:  redC,
		B:  redB,
		Tr: Tr,
		Pr: Pr,
	})

	return Tr * (delta*BH + delta*delta*CH), nil

}

// ResidualEntropy calculates the dimensionless residual entropy S^R/R using
// the second- and third-virial correlations.
//
// The result is evaluated as:
//
//	S^R/R = ln(Pr/(Delta*Tr)) - Delta*B_S - Delta^2*C_S/2
//
// It returns an error if Tr or Pr is not positive or if the reduced-density
// equation has no physically meaningful solution.
func ResidualEntropy(Tr, Pr, acentric float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}
	redB, err := abbott.ReducedB(Tr, acentric)
	if err != nil {
		return 0, err
	}
	redC, err := ReducedC(Tr, acentric)
	if err != nil {
		return 0, err
	}
	delta, err := Delta(DeltaArgs{
		C:  redC,
		B:  redB,
		Tr: Tr,
		Pr: Pr,
	})
	if err != nil {
		return 0, err
	}
	logTerm := math.Log(Pr / (Tr * delta))
	BS, err := bS(Tr, acentric)
	if err != nil {
		return 0, err
	}
	CS, err := cS(Tr, acentric)

	return logTerm - (delta*BS + delta*delta*CS/2), nil

}
