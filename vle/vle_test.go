package vle_test

import (
	"errors"
	"fmt"
	"strings"
	"testing"

	"github.com/rickykimani/zfactor/vle"
)

// TestSolverOptionDefaults checks the promise the type documents, that a
// zero value is usable.
//
// Every VLE calculation accepts these options, and most callers leave
// them unset, so the substituted defaults are what actually governs the
// iterations.
func TestSolverOptionDefaults(t *testing.T) {
	var zero vle.SolverOptions

	if got := zero.Tol(); got <= 0 {
		t.Errorf("default tolerance = %g; want a positive value", got)
	}

	if got := zero.MaxIter(); got <= 0 {
		t.Errorf("default iteration limit = %d; want a positive value", got)
	}

	// A non-positive setting is treated as unset rather than honoured,
	// since neither would terminate.
	negative := vle.SolverOptions{Tolerance: -1, MaxIterations: -1}

	if negative.Tol() != zero.Tol() {
		t.Errorf("a negative tolerance gives %g; want the default %g", negative.Tol(), zero.Tol())
	}

	if negative.MaxIter() != zero.MaxIter() {
		t.Errorf("a negative iteration limit gives %d; want the default %d", negative.MaxIter(), zero.MaxIter())
	}

	// A positive setting is honoured.
	explicit := vle.SolverOptions{Tolerance: 1e-14, MaxIterations: 7}

	if explicit.Tol() != 1e-14 {
		t.Errorf("tolerance = %g; want the supplied 1e-14", explicit.Tol())
	}

	if explicit.MaxIter() != 7 {
		t.Errorf("iteration limit = %d; want the supplied 7", explicit.MaxIter())
	}
}

// TestPhaseStateString checks that every state renders a name, since the
// single-phase error embeds one in its message.
func TestPhaseStateString(t *testing.T) {
	testCases := []struct {
		state vle.PhaseState
		want  string
	}{
		{vle.TwoPhase, "two-phase"},
		{vle.SubcooledLiquid, "subcooled liquid"},
		{vle.SuperheatedVapor, "superheated vapor"},
		{vle.PhaseState(99), "unknown"},
	}

	for _, tc := range testCases {
		if got := tc.state.String(); got != tc.want {
			t.Errorf("state %d renders as %q; want %q", int(tc.state), got, tc.want)
		}
	}
}

// TestSinglePhaseError checks that the error names the state and that the
// state survives wrapping.
//
// Distinguishing a subcooled liquid from a superheated vapor is the
// reason the type carries a state at all: a caller adjusting conditions
// to reach the two-phase region needs to know which way to move.
func TestSinglePhaseError(t *testing.T) {
	for _, state := range []vle.PhaseState{vle.SubcooledLiquid, vle.SuperheatedVapor} {
		err := &vle.SinglePhaseError{State: state}

		if got := err.Error(); !strings.Contains(got, state.String()) {
			t.Errorf("message %q does not name the state %q", got, state)
		}

		wrapped := fmt.Errorf("flashing the feed: %w", err)

		var single *vle.SinglePhaseError
		if !errors.As(wrapped, &single) {
			t.Fatalf("errors.As did not recover the error from %v", wrapped)
		}

		if single.State != state {
			t.Errorf("recovered state %v; want %v", single.State, state)
		}
	}
}
