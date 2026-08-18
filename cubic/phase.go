package cubic

import (
	"errors"
	"fmt"
	"math"
)

// Phase names which root of the cubic a property is being asked for.
//
// A state admits either one real root or three. Where it admits three, the
// smallest is the liquid and the largest the vapour, and the root between them
// is not a physical state. Where it admits one, every phase names that root:
// there is only one way for the fluid to be.
type Phase int

const (
	// StablePhase selects whichever phase is thermodynamically stable, which
	// is the one of lower fugacity. It is the zero value, so a caller who has
	// not thought about phase gets the physical answer rather than an
	// arbitrary one.
	//
	// Stability here means stable according to this equation of state: two
	// equations disagree about which phase exists near saturation, and both
	// disagree with the fluid. Near the critical point the two roots converge
	// and the comparison stops being meaningful before it stops being
	// computable.
	StablePhase Phase = iota

	// VaporPhase selects the largest real root, whether or not it is stable.
	VaporPhase

	// LiquidPhase selects the smallest real root, whether or not it is
	// stable. A metastable root is a legitimate thing to want a property
	// for, which is why naming a phase outright stays available.
	LiquidPhase
)

func (p Phase) String() string {
	switch p {
	case StablePhase:
		return "stable"
	case VaporPhase:
		return "vapor"
	case LiquidPhase:
		return "liquid"
	default:
		return fmt.Sprintf("Phase(%d)", int(p))
	}
}

// ErrNoRealRoot is returned when the equation of state has no real root at a
// state, so there is no volume for a property to be evaluated at.
var ErrNoRealRoot = errors.New("equation of state has no real root at this state")

// PhaseState is a solved state in the dimensionless form the property
// functions take.
//
// The three quantities always travel together — LogFugacity, ResidualEnthalpy
// and ResidualEntropy each take all of them — and forming them by hand means
// repeating the same three reductions at every call site.
type PhaseState struct {
	// V is the molar volume of the chosen phase, in the units implied by the
	// gas constant. It is the root itself, before any reduction.
	V float64

	// Z is the compressibility factor PV/(RT) of the chosen phase.
	Z float64

	// A is the dimensionless a P / (RT)².
	A float64

	// B is the dimensionless b P / (RT), written beta in the literature.
	B float64
}

// SolvePhase solves the equation of state at the temperature and pressure in
// cfg and returns the requested phase in dimensionless form.
//
// It is the step between SolveForVolume, which gives molar volumes, and the
// property functions, which take Z along with A and B.
func SolvePhase(cfg *EOSCfg, phase Phase) (*PhaseState, error) {
	if cfg == nil {
		return nil, errors.New("configuration error: config cannot be nil")
	}
	if cfg.T <= 0 {
		return nil, errors.New("configuration error: temperature must be greater than 0")
	}
	if cfg.P <= 0 {
		return nil, errors.New("configuration error: pressure must be greater than 0")
	}
	if cfg.R <= 0 {
		return nil, errors.New("configuration error: gas constant must be greater than 0")
	}

	result, err := SolveForVolume(cfg)
	if err != nil {
		return nil, err
	}

	roots := result.Clean()
	if len(roots) == 0 {
		return nil, fmt.Errorf("%w: T=%g P=%g", ErrNoRealRoot, cfg.T, cfg.P)
	}

	RT := cfg.R * cfg.T
	a := result.A * cfg.P / (RT * RT)
	b := result.B * cfg.P / RT

	// Clean returns the real roots in ascending order, so the ends are the
	// two phases. With a single root every phase names it.
	liquid, vapor := roots[0], roots[len(roots)-1]

	var volume float64

	switch phase {
	case LiquidPhase:
		volume = liquid
	case VaporPhase:
		volume = vapor
	case StablePhase:
		volume = stableRoot(cfg, liquid, vapor, a, b, RT)
	default:
		return nil, fmt.Errorf("unknown phase: %v", phase)
	}

	return &PhaseState{
		V: volume,
		Z: cfg.P * volume / RT,
		A: a,
		B: b,
	}, nil
}

// stableRoot returns whichever of the two outer roots is thermodynamically
// stable, which is the one of lower fugacity.
//
// This is equivalent to comparing the pressure against the saturation
// pressure, but does not need it: SaturationPressure is an iterative solve
// that can fail to converge, whereas the fugacity of each root follows
// directly from numbers already in hand. At the saturation pressure the two
// are equal by construction — that is the condition SaturationPressure exists
// to find — so the tie goes to the vapour by convention rather than to
// whichever way the comparison happens to fall.
func stableRoot(cfg *EOSCfg, liquid, vapor, a, b, RT float64) float64 {
	if liquid == vapor {
		return vapor
	}

	// A root at or below b leaves ln(Z - B) undefined, so it is no candidate.
	// This cannot happen for a root of the equation, but a root recovered
	// through floating point can land there.
	zLiquid := cfg.P * liquid / RT
	zVapor := cfg.P * vapor / RT

	switch {
	case zLiquid <= b:
		// Includes the case where neither is usable, which leaves nothing to
		// compare and no better answer than the vapour.
		return vapor
	case zVapor <= b:
		return liquid
	}

	lnPhiLiquid := LogFugacity(cfg, zLiquid, a, b)
	lnPhiVapor := LogFugacity(cfg, zVapor, a, b)

	// At the saturation pressure the two are equal in principle, but only to
	// the tolerance the pressure was found to, so the comparison would
	// otherwise be settled by the last bits of a converged iteration. Calling
	// anything inside that tolerance a tie makes the convention real.
	//
	// The band this covers is vanishingly narrow in pressure — the fugacities
	// separate by a few hundredths per bar for a wide two-phase region, so
	// equality to 1e-8 means within about a microbar of saturation — and so it
	// does not reach any state that is genuinely one phase.
	if math.Abs(lnPhiLiquid-lnPhiVapor) < fugacityTie {
		return vapor
	}

	if lnPhiLiquid < lnPhiVapor {
		return liquid
	}

	return vapor
}

// fugacityTie is the difference in ln phi below which the two phases are taken
// to be in equilibrium. It matches the convergence SaturationPressure works
// to, since that is what sets how close to equal the two can be found to be.
const fugacityTie = 1e-8
