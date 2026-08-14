package vle

import "fmt"

// PhaseState describes where a mixture sits relative to its two-phase
// region, which determines whether a flash calculation has an answer.
type PhaseState int

const (
	// TwoPhase means vapor and liquid coexist, so a flash has a solution.
	TwoPhase PhaseState = iota

	// SubcooledLiquid means the pressure is above the bubble pressure of
	// the feed, so nothing has evaporated.
	SubcooledLiquid

	// SuperheatedVapor means the pressure is below the dew pressure of
	// the feed, so nothing has condensed.
	SuperheatedVapor
)

func (p PhaseState) String() string {
	switch p {
	case TwoPhase:
		return "two-phase"
	case SubcooledLiquid:
		return "subcooled liquid"
	case SuperheatedVapor:
		return "superheated vapor"
	default:
		return "unknown"
	}
}

// SinglePhaseError reports that a mixture is wholly liquid or wholly
// vapor at the specified conditions, so no phase split exists to
// calculate.
//
// It describes a state rather than a failure: a feed separates only
// between its bubble and dew pressures, and asking for a flash outside
// that range is a reasonable question with a definite negative answer.
// The state is carried so that a caller can tell which side of the
// two-phase region the feed lies on without repeating the test:
//
//	var single *vle.SinglePhaseError
//	if errors.As(err, &single) {
//	    switch single.State {
//	    case vle.SubcooledLiquid: ...
//	    case vle.SuperheatedVapor: ...
//	    }
//	}
type SinglePhaseError struct {
	State PhaseState
}

func (e *SinglePhaseError) Error() string {
	return fmt.Sprintf(
		"the mixture is a %s at these conditions, so no phase split exists",
		e.State,
	)
}
