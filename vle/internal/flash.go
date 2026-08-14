package internal

import (
	"errors"
	"fmt"

	"github.com/rickykimani/zfactor/vle"
)

// ClassifyFeed reports where a feed of overall composition z sits
// relative to its two-phase region, given the equilibrium ratios K.
//
// The Rachford-Rice function evaluated at the ends of the vapor-fraction
// range supplies the answer without a separate bubble- or dew-point
// calculation:
//
//	F(0) = Σ zi Ki - 1
//	F(1) = 1 - Σ zi/Ki
//
// The first is positive only below the bubble pressure and the second
// negative only above the dew pressure, so both conditions together
// bracket the two-phase region.
func ClassifyFeed(z, k []float64) (vle.PhaseState, error) {
	if err := validateFlashInput(z, k); err != nil {
		return vle.TwoPhase, err
	}

	var atZero, atOne float64
	for i := range z {
		atZero += z[i] * k[i]
		atOne += z[i] / k[i]
	}

	switch {
	case atZero <= 1:
		// Σ zi Ki ≤ 1: the feed has not reached its bubble point.
		return vle.SubcooledLiquid, nil
	case atOne <= 1:
		// Σ zi/Ki ≤ 1: the feed is past its dew point.
		return vle.SuperheatedVapor, nil
	default:
		return vle.TwoPhase, nil
	}
}

// VaporFraction solves the Rachford-Rice equation for the fraction of a
// feed that is vapor at equilibrium.
//
// Writing the equation as the difference of the two summation conditions
//
//	F(V) = Σ zi(Ki - 1) / (1 + V(Ki - 1)) = 0
//
// has two advantages over either condition alone. Its derivative,
//
//	dF/dV = -Σ zi(Ki - 1)² / [1 + V(Ki - 1)]²,
//
// is negative wherever it is defined, so F decreases monotonically and
// the root is unique. And it is free of the spurious roots the separate
// forms carry: summing the vapor compositions gives an identity at
// V = 1, and the liquid compositions at V = 0, so a solver applied to
// either would be drawn to an endpoint that is not a phase split.
//
// The poles of F lie at V = 1/(1 - Ki), which for positive Ki is
// negative when Ki exceeds one and greater than one otherwise. None
// falls inside [0, 1], so F is smooth across the whole search interval
// and bisection cannot step onto a singularity.
//
// It returns a *vle.SinglePhaseError if the feed does not split at
// these conditions.
func VaporFraction(z, k []float64, opts vle.SolverOptions) (float64, error) {
	state, err := ClassifyFeed(z, k)
	if err != nil {
		return 0, err
	}

	if state != vle.TwoPhase {
		return 0, &vle.SinglePhaseError{State: state}
	}

	residual := func(v float64) (float64, error) {
		var sum float64
		for i := range z {
			sum += z[i] * (k[i] - 1) / (1 + v*(k[i]-1))
		}
		return sum, nil
	}

	// The classification above established that the residual changes sign
	// across the interval, so the bracket is guaranteed.
	return Bisection(residual, 0, 1, opts)
}

// FlashCompositions returns the equilibrium vapor and liquid
// compositions of a feed z split at vapor fraction v with equilibrium
// ratios k.
//
//	yi = zi Ki / (1 + V(Ki - 1))      xi = yi / Ki
//
// Both follow from the material balance zi = xi(1 - V) + yi V together
// with the definition Ki = yi/xi, so they satisfy that balance by
// construction.
func FlashCompositions(z, k []float64, v float64) (x, y []float64, err error) {
	if err := validateFlashInput(z, k); err != nil {
		return nil, nil, err
	}

	if v < 0 || v > 1 {
		return nil, nil, fmt.Errorf("vapor fraction %g lies outside [0, 1]", v)
	}

	x = make([]float64, len(z))
	y = make([]float64, len(z))

	for i := range z {
		y[i] = z[i] * k[i] / (1 + v*(k[i]-1))
		x[i] = y[i] / k[i]
	}

	return x, y, nil
}

// validateFlashInput checks the feed composition and equilibrium ratios.
//
// A non-positive equilibrium ratio would place a pole of the
// Rachford-Rice function inside the search interval, and is in any case
// not physical: it would mean a component present in one phase and
// wholly absent from the other.
func validateFlashInput(z, k []float64) error {
	if len(z) == 0 {
		return errors.New("no components provided")
	}

	if len(k) != len(z) {
		return fmt.Errorf(
			"got %d equilibrium ratios for %d components", len(k), len(z),
		)
	}

	var total float64

	for i, zi := range z {
		if zi < 0 || zi > 1 {
			return fmt.Errorf("component %d has an overall mole fraction of %g", i, zi)
		}

		if k[i] <= 0 {
			return fmt.Errorf("component %d has a non-positive equilibrium ratio of %g", i, k[i])
		}

		total += zi
	}

	if diff := total - 1; diff > flashTolerance || diff < -flashTolerance {
		return fmt.Errorf("overall mole fractions sum to %g; want 1", total)
	}

	return nil
}

// flashTolerance is the departure from unity allowed of the overall
// composition, matching the tolerance the VLE packages use elsewhere.
const flashTolerance = 1e-6
