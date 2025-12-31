// Package abbott provides the generalized correlations for the second virial coefficient
//
// The second virial coefficient B is calculated as:
//
//	B * Pc / (R * Tc) = B0 + ω * B1 
//
// where B0 is the simple fluid contribution and B1 is the correction term for
// acentric factor (ω).
package abbott

import (
	"math"

	"github.com/rickykimani/zfactor"
)

// B0 calculates the simple fluid contribution to the second virial coefficient.
//
//	B0 = 0.083 - 0.422 / Tr^1.6
//
// It returns an error if Tr <= 0.
func B0(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return 0.083 - 0.422/math.Pow(Tr, 1.6), nil
}

// B1 calculates the correction term for the second virial coefficient based on the acentric factor.
//
//	B1 = 0.139 - 0.172 / Tr^4.2
//
// It returns an error if Tr <= 0.
func B1(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return 0.139 - 0.172/math.Pow(Tr, 4.2), nil
}

// DB0 calculates the first derivative of B0 with respect to reduced temperature (Tr).
//
//	dB0/dTr = 0.675 / Tr^2.6
//
// It returns an error if Tr <= 0.
func DB0(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return 0.675 / math.Pow(Tr, 2.6), nil
}

// DB1 calculates the first derivative of B1 with respect to reduced temperature (Tr).
//
//	dB1/dTr = 0.722 / Tr^5.2
//
// It returns an error if Tr <= 0.
func DB1(Tr float64) (float64, error) {
	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	return 0.722 / math.Pow(Tr, 5.2), nil
}
