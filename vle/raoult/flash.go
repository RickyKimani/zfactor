package raoult

import (
	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/vle"
	"github.com/rickykimani/zfactor/vle/internal"
)

// FlashInput supplies the information required for an isothermal,
// isobaric flash calculation.
//
// The composition is the overall composition of the feed, which is
// neither the liquid nor the vapor composition but the average of the two
// weighted by their amounts.
type FlashInput interface {
	Composition() []float64
	Pressure() float64
	PSat() ([]float64, error)
	SolverOptions() vle.SolverOptions
}

// FlashResult contains the amounts and compositions of the two phases a
// feed separates into.
//
// The amounts are per mole of feed, so they sum to one.
type FlashResult struct {
	V float64   // moles of vapor per mole of feed
	L float64   // moles of liquid per mole of feed
	X []float64 // equilibrium liquid composition
	Y []float64 // equilibrium vapor composition
}

// FlashPT calculates the phase split of a feed at a specified
// temperature and pressure using Raoult's law.
//
// A liquid at or above its bubble pressure partially evaporates when the
// pressure is reduced below it, and the calculation determines how much
// of the feed ends up in each phase and what each phase contains. Duhem's
// theorem fixes the state: for a feed of known composition, specifying
// two variables determines the rest.
//
// Under Raoult's law the equilibrium ratios
//
//	Ki = yi/xi = Pi_sat / P
//
// depend only on the temperature and pressure, so they are known before
// the split is and a single solution of the Rachford-Rice equation
// suffices. Where activity coefficients are involved the ratios depend on
// the liquid composition being solved for, and the calculation must
// iterate; see the modified Raoult's law package.
//
// A feed only separates between its bubble and dew pressures. Outside
// that range the result is a *vle.SinglePhaseError naming which side the
// feed falls on, which describes the state rather than reporting a
// failure of the calculation.
func FlashPT(input FlashInput) (FlashResult, error) {
	// The feed composition and saturation pressures are validated the
	// same way as for a bubble- or dew-point calculation.
	res, err := preparePressureInput(input)
	if err != nil {
		return FlashResult{}, err
	}

	p := input.Pressure()
	if p <= 0 {
		return FlashResult{}, zfactor.ErrPressure
	}

	z := res.comp

	k := make([]float64, res.n)
	for i, psat := range res.psat {
		k[i] = psat / p
	}

	opts := input.SolverOptions()

	v, err := internal.VaporFraction(z, k, opts)
	if err != nil {
		return FlashResult{}, err
	}

	x, y, err := internal.FlashCompositions(z, k, v)
	if err != nil {
		return FlashResult{}, err
	}

	return FlashResult{
		V: v,
		L: 1 - v,
		X: x,
		Y: y,
	}, nil
}
