package modified_raoult

import (
	"errors"
	"fmt"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/vle"
	"github.com/rickykimani/zfactor/vle/internal"
)

// ErrNoAzeotrope reports that the specified binary system does not form
// an azeotrope under the given conditions.
//
// This is a normal result rather than a computational failure: many
// binary mixtures simply exhibit no composition at which the relative
// volatility reaches unity. Callers should distinguish it with
// errors.Is.
var ErrNoAzeotrope = errors.New(
	"no azeotrope exists for the specified binary system",
)

const (
	// edge is the offset from the pure-component limits at which the
	// composition sweep begins and ends. Activity models are evaluated
	// at x = edge and x = 1-edge rather than exactly at 0 and 1, where
	// some correlations are ill-conditioned.
	edge = 1e-6

	// scanIntervals is the number of subintervals into which the
	// composition range is divided when searching for azeotropes.
	//
	// A binary mixture may form more than one azeotrope, so the search
	// cannot rely on the residual changing sign between the two
	// pure-component limits: with an even number of roots it does not.
	// The range is therefore swept and every subinterval containing a
	// sign change is refined separately.
	//
	// The sweep resolves azeotropes separated by more than roughly
	// 1/scanIntervals in mole fraction. Raising it improves resolution
	// at proportional cost, which is most noticeable for AzeotropeT
	// where every evaluation solves a bubble-temperature problem.
	scanIntervals = 50
)

// AzeotropeResult contains one azeotropic state of a binary mixture.
//
// At an azeotrope the liquid and vapor compositions are identical, so X
// and Y hold the same values and are both reported for symmetry with
// the bubble- and dew-point results.
type AzeotropeResult struct {
	X []float64 // azeotropic liquid composition
	Y []float64 // equilibrium vapor composition (equal to X)
	T float64   // azeotropic temperature (°C)
	P float64   // azeotropic pressure (kPa)
}

// AzeotropePInput supplies the information required for an isothermal
// azeotrope calculation.
//
// Unlike the bubble- and dew-point inputs, no composition is supplied:
// the azeotropic composition is the quantity being solved for. The
// input types of this package satisfy this interface, so the
// composition field may simply be left unset.
type AzeotropePInput interface {
	Temperature() float64
	PSat() ([]float64, error)
	ActivityModel() ActivityModel
	SolverOptions() vle.SolverOptions
}

// AzeotropeTInput supplies the information required for an isobaric
// azeotrope calculation.
//
// Saturation pressures are supplied as Antoine correlations rather than
// fixed values because the azeotropic temperature is unknown and varies
// throughout the search.
type AzeotropeTInput interface {
	Pressure() float64
	AntoineModels() []antoine.Model
	ActivityModel() ActivityModel
	SolverOptions() vle.SolverOptions
}

// AzeotropeP calculates the azeotropic compositions and pressures of a
// binary mixture at a specified temperature using modified Raoult's law.
//
// At an azeotrope the liquid and vapor compositions coincide (x = y),
// which reduces the equilibrium relation yi P = xi γi Pi_sat to
//
//	P = γi Pi_sat   for every component,
//
// equivalently a relative volatility of unity:
//
//	α₁₂ = γ₁ P₁_sat / (γ₂ P₂_sat) = 1.
//
// Each azeotropic composition is therefore a root of
//
//	f(x₁) = γ₁ P₁_sat - γ₂ P₂_sat
//
// on the open interval (0, 1). The difference form is used in place of
// the volatility ratio to avoid dividing by a vanishing activity
// coefficient or saturation pressure.
//
// A binary mixture may form more than one azeotrope, so every root is
// reported, ordered by increasing x₁. If the mixture forms none,
// ErrNoAzeotrope is returned.
//
// Only binary mixtures are supported.
func AzeotropeP(input AzeotropePInput) ([]AzeotropeResult, error) {

	T := input.Temperature()
	if T <= -tempConv {
		return nil, zfactor.ErrTemp
	}

	psat, err := input.PSat()
	if err != nil {
		return nil, err
	}

	if err := validateBinary(len(psat), "saturation pressures"); err != nil {
		return nil, err
	}

	if err := validatePSat(psat); err != nil {
		return nil, err
	}

	activityModel := input.ActivityModel()
	if activityModel == nil {
		return nil, errors.New("activity model is nil")
	}

	// Residual whose roots are the azeotropic compositions.
	residual := func(x1 float64) (float64, error) {
		gamma, err := activityCoefficients(activityModel, T, binary(x1))
		if err != nil {
			return 0, err
		}
		return gamma[0]*psat[0] - gamma[1]*psat[1], nil
	}

	compositions, err := azeotropicCompositions(residual, input.SolverOptions())
	if err != nil {
		return nil, err
	}

	results := make([]AzeotropeResult, 0, len(compositions))

	for _, x1 := range compositions {
		x := binary(x1)

		gamma, err := activityCoefficients(activityModel, T, x)
		if err != nil {
			return nil, err
		}

		results = append(results, AzeotropeResult{
			X: x,
			Y: binary(x1),
			T: T,
			P: gamma[0] * psat[0],
		})
	}

	return results, nil
}

// AzeotropeT calculates the azeotropic compositions and temperatures of
// a binary mixture at a specified pressure using modified Raoult's law.
//
// The azeotropic condition is the same as for AzeotropeP, but the
// temperature is no longer known in advance: both the activity
// coefficients and the saturation pressures shift as it changes. The
// calculation is therefore nested. For each trial composition the
// bubble temperature at the specified pressure is obtained from
// BubbleT, and the azeotropic residual
//
//	f(x₁) = γ₁ P₁_sat(T) - γ₂ P₂_sat(T)
//
// is evaluated at that temperature. The outer search locates the
// compositions at which the residual vanishes, which are the states
// satisfying both the bubble-point and azeotropic conditions.
//
// This reuses the tested bubble-temperature solver at the cost of an
// inner iteration per residual evaluation; a simultaneous solution of
// both conditions would be faster but considerably less transparent.
//
// Every azeotrope is reported, ordered by increasing x₁. If the mixture
// forms none, ErrNoAzeotrope is returned.
//
// Only binary mixtures are supported.
func AzeotropeT(input AzeotropeTInput) ([]AzeotropeResult, error) {

	P := input.Pressure()
	if P <= 0 {
		return nil, zfactor.ErrPressure
	}

	models := input.AntoineModels()

	if err := validateBinary(len(models), "Antoine models"); err != nil {
		return nil, err
	}

	for i, model := range models {
		if model == nil {
			return nil, fmt.Errorf("antoine model %d is nil", i)
		}
	}

	activityModel := input.ActivityModel()
	if activityModel == nil {
		return nil, errors.New("activity model is nil")
	}

	opts := input.SolverOptions()

	// bubbleTemperature returns the bubble temperature (°C) of the
	// specified liquid composition at the specified pressure.
	bubbleTemperature := func(x []float64) (float64, error) {
		res, err := BubbleT(MixtureInput{
			P:            P,
			Compositions: x,
			Antoine:      models,
			Activity:     activityModel,
			Options:      opts,
		})
		if err != nil {
			return 0, err
		}
		return res.T, nil
	}

	// Residual evaluated along the bubble-point locus.
	residual := func(x1 float64) (float64, error) {
		x := binary(x1)

		T, err := bubbleTemperature(x)
		if err != nil {
			return 0, err
		}

		gamma, err := activityCoefficients(activityModel, T, x)
		if err != nil {
			return 0, err
		}

		psat := make([]float64, len(models))
		for i, model := range models {
			psat[i], err = saturationPressure(model, T)
			if err != nil {
				return 0, err
			}
		}

		return gamma[0]*psat[0] - gamma[1]*psat[1], nil
	}

	compositions, err := azeotropicCompositions(residual, opts)
	if err != nil {
		return nil, err
	}

	results := make([]AzeotropeResult, 0, len(compositions))

	for _, x1 := range compositions {
		x := binary(x1)

		T, err := bubbleTemperature(x)
		if err != nil {
			return nil, err
		}

		results = append(results, AzeotropeResult{
			X: x,
			Y: binary(x1),
			T: T,
			P: P,
		})
	}

	return results, nil
}

// azeotropicCompositions locates every composition at which the
// supplied azeotropic residual vanishes.
//
// The composition range is swept in scanIntervals steps and each
// subinterval across which the residual changes sign is refined by
// bisection, so the iteration cannot leave the physically meaningful
// interval (0, 1).
//
// Sweeping is necessary because a binary mixture may form more than one
// azeotrope. Comparing only the two pure-component limits detects a
// root solely when their signs differ, which fails for any even number
// of roots: a double azeotrope would be reported as no azeotrope at all.
//
// Roots are returned in increasing order. ErrNoAzeotrope is returned
// when the sweep finds none.
func azeotropicCompositions(
	residual func(float64) (float64, error),
	opts vle.SolverOptions,
) ([]float64, error) {

	lo, hi := edge, 1-edge
	step := (hi - lo) / scanIntervals

	var roots []float64

	// Two roots closer together than the solver tolerance are
	// indistinguishable, so a root already recorded is not repeated.
	record := func(x float64) {
		for _, r := range roots {
			if math.Abs(r-x) <= opts.Tol() {
				return
			}
		}
		roots = append(roots, x)
	}

	prevX := lo
	prevF, err := residual(prevX)
	if err != nil {
		return nil, err
	}
	if prevF == 0 {
		record(prevX)
	}

	for i := 1; i <= scanIntervals; i++ {
		x := lo + step*float64(i)
		if i == scanIntervals {
			x = hi
		}

		f, err := residual(x)
		if err != nil {
			return nil, err
		}

		switch {
		case f == 0:
			// Root sits exactly on a sweep point.
			record(x)

		case prevF != 0 && math.Signbit(f) != math.Signbit(prevF):
			root, err := internal.Bisection(residual, prevX, x, opts)
			if err != nil {
				return nil, err
			}
			record(root)
		}

		prevX, prevF = x, f
	}

	if len(roots) == 0 {
		return nil, ErrNoAzeotrope
	}

	return roots, nil
}

// binary returns the two-component composition vector [x1, 1-x1].
func binary(x1 float64) []float64 {
	return []float64{x1, 1 - x1}
}

// validateBinary reports whether a supplied slice describes exactly two
// components, naming the quantity in the error message.
func validateBinary(n int, quantity string) error {
	if n != 2 {
		return fmt.Errorf(
			"azeotrope calculations support binary mixtures only: got %d %s",
			n, quantity,
		)
	}
	return nil
}
