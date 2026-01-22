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
