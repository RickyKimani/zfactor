package ov

import (
	"errors"
	"math"

	"github.com/rickykimani/zfactor/abbott"
)

// bH combines the Abbott second-virial correlation and its derivative for the
// residual-enthalpy expression:
//
//	B_H = (B0 - Tr*dB0/dTr) + omega*(B1 - Tr*dB1/dTr)
func bH(Tr, acentric float64) (float64, error) {
	B0, err := abbott.B0(Tr)
	if err != nil {
		return 0, err
	}
	B1, err := abbott.B1(Tr)
	if err != nil {
		return 0, err
	}
	DB0, err := abbott.DB0(Tr)
	if err != nil {
		return 0, err
	}
	DB1, err := abbott.DB1(Tr)
	if err != nil {
		return 0, err
	}

	return (B0 - Tr*DB0) + acentric*(B1-Tr*DB1), nil

}

// cH combines the Orbey-Vera third-virial correlation and its derivative for
// the residual-enthalpy expression. The third-virial derivative carries the
// one-half factor from the density-integral formulation:
//
//	C_H = (C0 - Tr*dC0/(2*dTr)) + omega*(C1 - Tr*dC1/(2*dTr))
func cH(Tr, acentric float64) (float64, error) {
	C0, err := C0(Tr)
	if err != nil {
		return 0, err
	}
	C1, err := C1(Tr)
	if err != nil {
		return 0, err
	}
	DC0, err := DC0(Tr)
	if err != nil {
		return 0, err
	}
	DC1, err := DC1(Tr)
	if err != nil {
		return 0, err
	}

	return (C0 - Tr*DC0/2) + acentric*(C1-Tr*DC1/2), nil

}

// bS combines the Abbott second-virial correlation and its derivative for the
// residual-entropy expression:
//
//	B_S = (B0 + Tr*dB0/dTr) + omega*(B1 + Tr*dB1/dTr)
func bS(Tr, acentric float64) (float64, error) {
	B0, err := abbott.B0(Tr)
	if err != nil {
		return 0, err
	}
	B1, err := abbott.B1(Tr)
	if err != nil {
		return 0, err
	}
	DB0, err := abbott.DB0(Tr)
	if err != nil {
		return 0, err
	}
	DB1, err := abbott.DB1(Tr)
	if err != nil {
		return 0, err
	}

	return (B0 + Tr*DB0) + acentric*(B1+Tr*DB1), nil

}

// cS combines the Orbey-Vera third-virial correlation and its derivative for
// the residual-entropy expression:
//
//	C_S = (C0 + Tr*dC0/dTr) + omega*(C1 + Tr*dC1/dTr)
func cS(Tr, acentric float64) (float64, error) {
	C0, err := C0(Tr)
	if err != nil {
		return 0, err
	}
	C1, err := C1(Tr)
	if err != nil {
		return 0, err
	}
	DC0, err := DC0(Tr)
	if err != nil {
		return 0, err
	}
	DC1, err := DC1(Tr)
	if err != nil {
		return 0, err
	}

	return (C0 + Tr*DC0) + acentric*(C1+Tr*DC1), nil

}

// clean returns the smallest positive real root from the delta equation.
// Complex roots and non-positive roots are ignored. An error is returned if
// no physically meaningful root is found.
func clean(deltas [3]complex128) (float64, error) {
	const imagTol = 1e-9

	var (
		root  float64
		found bool
	)

	for _, delta := range deltas {
		if math.Abs(imag(delta)) >= imagTol {
			continue
		}

		r := real(delta)
		if r <= 0 {
			continue
		}

		if !found || r < root {
			root = r
			found = true
		}
	}

	if !found {
		return 0, errors.New("no physically meaningful positive real roots found")
	}

	return root, nil
}
