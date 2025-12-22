package virial

import (
	"github.com/rickykimani/zfactor"
)

// SolveForVolumeTwoTerm solves the 2-term virial equation for molar volume.
// It uses the approximation V = RT/P + B.
func SolveForVolumeTwoTerm(T, P, R, B float64) (float64, error) {
	if P <= 0 {
		return 0, zfactor.ErrPressure
	}
	if P > 15 {
		return 0, zfactor.ErrHighPressureTwoTerm
	}
	if T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if R <= 0 {
		return 0, zfactor.ErrUniversalConst
	}
	if B == 0 {
		return 0, zfactor.ErrVirialCoeff
	}

	return (R * T / P) + B, nil
}

// SolveForVolumeThreeTerm solves the 3-term virial equation (Leiden form) for molar volume.
// The equation is Z = 1 + B/V + C/V^2, which rearranges to a cubic equation in V.
func SolveForVolumeThreeTerm(T, P, R, B, C float64) ([3]complex128, error) {
	if P <= 0 {
		return [3]complex128{}, zfactor.ErrPressure
	}
	if T <= 0 {
		return [3]complex128{}, zfactor.ErrTemp
	}
	if R <= 0 {
		return [3]complex128{}, zfactor.ErrUniversalConst
	}
	if B == 0 || C == 0 {
		return [3]complex128{}, zfactor.ErrVirialCoeff
	}

	a := P / (R * T)
	b := -1.0
	c := -B
	d := -C

	return zfactor.SolveCubic(a, b, c, d)
}

// CompressibilityTwoTerm calculates the compressibility factor Z using the 2-term virial equation.
// Z = 1 + BP/RT
func CompressibilityTwoTerm(T, P, R, B float64) (float64, error) {
	if P <= 0 {
		return 0, zfactor.ErrPressure
	}
	if P > 15 {
		return 0, zfactor.ErrHighPressureTwoTerm
	}
	if T <= 0 {
		return 0, zfactor.ErrTemp
	}
	if R <= 0 {
		return 0, zfactor.ErrUniversalConst
	}
	if B == 0 {
		return 0, zfactor.ErrVirialCoeff
	}

	return 1 + (B*P)/(R*T), nil
}

// CompressibilityThreeTerm calculates the compressibility factor Z using the 3-term virial equation.
// Z = 1 + B/V + C/V^2
func CompressibilityThreeTerm(V, B, C float64) (float64, error) {
	if V <= 0 {
		return 0, zfactor.ErrVolume
	}
	if B == 0 || C == 0 {
		return 0, zfactor.ErrVirialCoeff
	}

	return 1 + B/V + C/(V*V), nil
}
