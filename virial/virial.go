// Package virial implements the virial equations of state, which express
// the compressibility factor as a power series in density.
//
// The series has a basis in statistical mechanics rather than being an
// empirical fit: the second coefficient B accounts for interactions
// between pairs of molecules, the third coefficient C for interactions
// among triples, and so on. Truncating it after a given term is
// therefore a statement about how dense the fluid is, and the equations
// here are for gases rather than liquids.
//
// Two truncations are provided. The two-term form
//
//	Z = 1 + BP/(RT)
//
// is linear in pressure and is the one obtained by truncating the
// pressure series after B. It is reliable only at low pressure, and the
// functions using it refuse to evaluate above 15 bar rather than return
// a number the equation does not support.
//
// The three-term form is written in the density, or Leiden, series
//
//	Z = 1 + B/V + C/V^2,
//
// which rearranges to a cubic in the molar volume and so is solved with
// the same root finder as the cubic equations of state. Retaining the
// third coefficient extends the usable range to moderate densities.
//
// Coefficients must be supplied by the caller through zfactor.Args.
// Generalized correlations for B are available in the abbott package,
// which estimates it from the reduced temperature and acentric factor.
//
// Units follow whatever is used consistently for the gas constant: with
// R in bar·cm³/(mol·K), pressures are in bar, temperatures in kelvin,
// B in cm³/mol and C in cm⁶/mol².
package virial

import (
	"github.com/rickykimani/zfactor"
)

// SolveForVolumeTwoTerm solves the 2-term virial equation for molar volume.
// It uses the approximation V = RT/P + B.
//
// Required Args:
//   - T: Temperature
//   - P: Pressure
//   - R: Gas Constant
//   - B: Second virial coefficient
func SolveForVolumeTwoTerm(args zfactor.Args) (float64, error) {
	if args.P <= 0 {
		return 0, zfactor.ErrPressure
	}
	if args.P > 15 {
		return 0, zfactor.ErrHighPressureTwoTerm
	}
	if args.T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if args.R <= 0 {
		return 0, zfactor.ErrUniversalConst
	}
	if args.B == 0 {
		return 0, zfactor.ErrVirialCoeff
	}

	return (args.R * args.T / args.P) + args.B, nil
}

// SolveForVolumeThreeTerm solves the 3-term virial equation (Leiden form) for molar volume.
// The equation is Z = 1 + B/V + C/V^2, which rearranges to a cubic equation in V.
//
// Required Args:
//   - T: Temperature
//   - P: Pressure
//   - R: Gas Constant
//   - B: Second virial coefficient
//   - C: Third virial coefficient
func SolveForVolumeThreeTerm(args zfactor.Args) ([3]complex128, error) {
	if args.P <= 0 {
		return [3]complex128{}, zfactor.ErrPressure
	}
	if args.T <= 0 {
		return [3]complex128{}, zfactor.ErrTemp
	}
	if args.R <= 0 {
		return [3]complex128{}, zfactor.ErrUniversalConst
	}
	if args.B == 0 || args.C == 0 {
		return [3]complex128{}, zfactor.ErrVirialCoeff
	}

	a := args.P / (args.R * args.T)
	b := -1.0
	c := -args.B
	d := -args.C

	return zfactor.SolveCubic(a, b, c, d)
}

// CompressibilityTwoTerm calculates the compressibility factor Z using the 2-term virial equation.
// Z = 1 + BP/RT
//
// Required Args:
//   - T: Temperature
//   - P: Pressure
//   - R: Gas Constant
//   - B: Second virial coefficient
func CompressibilityTwoTerm(args zfactor.Args) (float64, error) {
	if args.P <= 0 {
		return 0, zfactor.ErrPressure
	}
	if args.P > 15 {
		return 0, zfactor.ErrHighPressureTwoTerm
	}
	if args.T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if args.R <= 0 {
		return 0, zfactor.ErrUniversalConst
	}
	if args.B == 0 {
		return 0, zfactor.ErrVirialCoeff
	}

	return 1 + (args.B*args.P)/(args.R*args.T), nil
}

// CompressibilityThreeTerm calculates the compressibility factor Z using the 3-term virial equation.
// Z = 1 + B/V + C/V^2
//
// Required Args:
//   - B: Second virial coefficient
//   - C: Third virial coefficient
func CompressibilityThreeTerm(V float64, args zfactor.Args) (float64, error) {
	if V <= 0 {
		return 0, zfactor.ErrVolume
	}
	if args.B == 0 || args.C == 0 {
		return 0, zfactor.ErrVirialCoeff
	}

	return 1 + args.B/V + args.C/(V*V), nil
}
