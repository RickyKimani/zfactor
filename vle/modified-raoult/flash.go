package modified_raoult

import (
	"errors"

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
	Temperature() float64
	Pressure() float64
	Composition() []float64
	PSat() ([]float64, error)
	ActivityModel() ActivityModel
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
// temperature and pressure using the modified Raoult's law.
//
// The equilibrium ratios
//
//	Ki = yi/xi = γi Pi_sat / P
//
// depend on the activity coefficients, and those depend on the liquid
// composition the calculation is solving for. Unlike the ideal case the
// ratios are therefore not known in advance, and the calculation
// alternates between two steps: solving the Rachford-Rice equation for
// the vapor fraction at the current ratios, and re-evaluating the
// activity coefficients at the liquid composition that results. It
// converges when the compositions stop moving.
//
// A bubble- and a dew-pressure calculation at the feed composition
// precede the iteration. They establish the range of pressures over which
// the feed separates at all, and they supply the starting point: the
// vapor fraction is interpolated linearly between the two, and the
// activity coefficients along with it.
//
// Both are necessary rather than merely convenient. The bounds of the
// two-phase region cannot be read from the equilibrium ratios at the feed
// composition, because those ratios depend on the composition of the
// liquid, which near the dew point differs greatly from the feed. Under
// Raoult's law the ratios are independent of composition and the bounds
// do follow from the feed alone; that shortcut is not available here.
//
// Outside the range the result is a *vle.SinglePhaseError naming which
// side the feed falls on, which describes the state rather than reporting
// a failure of the calculation.
//
// The vapor phase is taken to be ideal, as it is throughout this
// package: the full formulation carries a factor Φi accounting for
// vapor-phase non-ideality and the Poynting correction, and here that
// factor is unity. The approximation is the one that distinguishes the
// modified Raoult's law from the general gamma/phi formulation, and it
// holds at low pressure.
func FlashPT(input FlashInput) (FlashResult, error) {
	t := input.Temperature()
	if t <= -tempConv {
		return FlashResult{}, zfactor.ErrTemp
	}

	p := input.Pressure()
	if p <= 0 {
		return FlashResult{}, zfactor.ErrPressure
	}

	z := input.Composition()
	if err := validateComposition(z); err != nil {
		return FlashResult{}, err
	}

	psat, err := input.PSat()
	if err != nil {
		return FlashResult{}, err
	}

	if err := validatePSat(psat); err != nil {
		return FlashResult{}, err
	}

	if len(psat) != len(z) {
		return FlashResult{}, errors.New(
			"number of saturation pressures must match number of components",
		)
	}

	activityModel := input.ActivityModel()
	if activityModel == nil {
		return FlashResult{}, errors.New("activity model is nil")
	}

	opts := input.SolverOptions()

	// ratios returns the equilibrium ratios for a given liquid
	// composition.
	ratios := func(x []float64) ([]float64, error) {
		gamma, err := activityCoefficients(activityModel, t, x)
		if err != nil {
			return nil, err
		}

		k := make([]float64, len(x))
		for i := range x {
			k[i] = gamma[i] * psat[i] / p
		}

		return k, nil
	}

	// Establish the two-phase range and the starting point.
	bubble, err := BubbleP(SaturationPressureInput{
		T: t, Compositions: z, PSats: psat, Activity: activityModel, Options: opts,
	})
	if err != nil {
		return FlashResult{}, err
	}

	dew, err := DewP(SaturationPressureInput{
		T: t, Compositions: z, PSats: psat, Activity: activityModel, Options: opts,
	})
	if err != nil {
		return FlashResult{}, err
	}

	switch {
	case p > bubble.P:
		return FlashResult{}, &vle.SinglePhaseError{State: vle.SubcooledLiquid}
	case p < dew.P:
		return FlashResult{}, &vle.SinglePhaseError{State: vle.SuperheatedVapor}
	}

	// Interpolate the liquid composition between the two boundaries. At
	// the bubble point the liquid is the feed; at the dew point it is the
	// liquid the dew calculation found.
	fraction := 0.0
	if span := bubble.P - dew.P; span > 0 {
		fraction = (bubble.P - p) / span
	}

	x := make([]float64, len(z))
	for i := range z {
		x[i] = z[i] + fraction*(dew.X[i]-z[i])
	}
	x = normalise(x)

	// Each pass solves for the vapor fraction at the current ratios and
	// returns the liquid composition it implies, which the next pass uses
	// to re-evaluate the activity coefficients. Successive substitution on
	// that composition is what converges.
	converged, err := internal.FixedPoint(x, func(current []float64) ([]float64, error) {
		k, err := ratios(current)
		if err != nil {
			return nil, err
		}

		v, err := internal.VaporFraction(z, k, opts)
		if err != nil {
			return nil, err
		}

		next, _, err := internal.FlashCompositions(z, k, v)
		if err != nil {
			return nil, err
		}

		// The compositions sum to one only at the exact root, and the
		// vapor fraction is known to the solver tolerance, so they are
		// renormalised before the activity model sees them. That model
		// validates its composition against a tighter tolerance than the
		// flash converges to, and would otherwise reject a composition
		// that is correct to within the accuracy asked of it.
		return normalise(next), nil
	}, opts)
	if err != nil {
		return FlashResult{}, err
	}

	// Re-solve at the converged composition so that the vapor fraction and
	// both compositions describe one state. The value carried out of the
	// final iteration was found from the composition before it.
	k, err := ratios(converged)
	if err != nil {
		return FlashResult{}, err
	}

	v, err := internal.VaporFraction(z, k, opts)
	if err != nil {
		return FlashResult{}, err
	}

	liquid, vapor, err := internal.FlashCompositions(z, k, v)
	if err != nil {
		return FlashResult{}, err
	}

	return FlashResult{
		V: v,
		L: 1 - v,
		X: liquid,
		Y: vapor,
	}, nil
}

// normalise scales a composition vector to sum to one.
//
// It is used between iterations rather than on the result: the material
// balance zi = xi(1 - V) + yi V holds exactly for any vapor fraction,
// whereas the compositions sum to one only when that fraction is the
// exact root, so scaling the returned values would trade an exact
// balance for a tidier sum.
func normalise(x []float64) []float64 {
	var total float64
	for _, xi := range x {
		total += xi
	}

	if total == 0 {
		return x
	}

	scaled := make([]float64, len(x))
	for i, xi := range x {
		scaled[i] = xi / total
	}

	return scaled
}
