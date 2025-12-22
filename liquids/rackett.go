package liquids

import (
	"math"

	"github.com/rickykimani/zfactor"
)

// Vsat calculates the saturated liquid molar volume using the Rackett equation.
//
// The Rackett equation is given by:
//
//	Vsat = Vc * Zc^((1 - Tr)^(2/7))
//
// Where:
//   - Vc is the critical molar volume.
//   - Zc is the critical compressibility factor.
//   - Tr is the reduced temperature (T/Tc).
//
// This correlation is typically accurate to within 1-2% for non-polar fluids.
func Vsat(Vc, Zc, Tr float64) (float64, error) {
	if Vc <= 0 || Zc <= 0 {
		return 0, zfactor.ErrCriticalProp
	}

	if Tr <= 0 {
		return 0, zfactor.ErrInvalidTr
	}

	// Big brain move: Square (1-Tr) first to avoid NaN when Tr > 1.
	v := Vc * math.Pow(Zc, math.Pow((1-Tr)*(1-Tr), 1.0/7.0))

	return v, nil
}
