package state_test

import (
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/state"
	"github.com/rickykimani/zfactor/substance"
)

// TestNewState checks the guards on the conditions a state is built
// from, and that a valid one carries them through.
func TestNewState(t *testing.T) {
	t.Run("valid", func(t *testing.T) {
		got, err := state.NewState(substance.Ethane, 299, 32)
		if err != nil {
			t.Fatalf("NewState returned an unexpected error: %v", err)
		}

		if got.Substance != substance.Ethane {
			t.Errorf("substance = %v; want ethane", got.Substance)
		}

		if got.Temperature != 299 || got.Pressure != 32 {
			t.Errorf("conditions = (%g, %g); want (299, 32)", got.Temperature, got.Pressure)
		}
	})

	invalid := []struct {
		name string
		T, P float64
	}{
		{"zero temperature", 0, 32},
		{"negative temperature", -10, 32},
		{"zero pressure", 299, 0},
		{"negative pressure", 299, -5},
	}

	for _, tc := range invalid {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := state.NewState(substance.Ethane, tc.T, tc.P); err == nil {
				t.Error("expected an error; got nil")
			}
		})
	}
}

// TestDrawPVRejectsMismatchedSubstances checks that a diagram may only be
// drawn for states of one substance.
//
// The axes are scaled from the substance's critical properties and the
// saturation dome drawn from its equation of state, so plotting two
// species on one figure would place them on incompatible axes.
func TestDrawPVRejectsMismatchedSubstances(t *testing.T) {
	ethane, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	propane, err := state.NewState(substance.Propane, 299, 32)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	cfg := &state.PVConfig{Type: &cubic.PR{}}
	output := filepath.Join(t.TempDir(), "mixed.svg")

	err = state.DrawPV(cfg, output, ethane, propane)
	if err == nil {
		t.Fatal("expected an error for states of different substances; got nil")
	}

	if !strings.Contains(err.Error(), "mismatch") {
		t.Errorf("error = %q; want it to name the mismatch", err)
	}
}

// TestDrawPVRequiresAState checks that a diagram cannot be drawn from no
// states, since the substance is taken from them.
func TestDrawPVRequiresAState(t *testing.T) {
	cfg := &state.PVConfig{Type: &cubic.PR{}}

	if err := state.DrawPV(cfg, filepath.Join(t.TempDir(), "empty.svg")); err == nil {
		t.Error("expected an error; got nil")
	}
}

// TestDrawPVConfigurationErrors checks the guards on the configuration
// itself.
func TestDrawPVConfigurationErrors(t *testing.T) {
	ethane, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	output := filepath.Join(t.TempDir(), "diagram.svg")

	t.Run("nil configuration", func(t *testing.T) {
		if err := state.DrawPV(nil, output, ethane); err == nil {
			t.Error("expected an error; got nil")
		}
	})

	t.Run("no equation of state", func(t *testing.T) {
		if err := state.DrawPV(&state.PVConfig{}, output, ethane); err == nil {
			t.Error("expected an error; got nil")
		}
	})
}

// TestDrawPVSuggestsTheNearestExtension checks the guidance offered when
// the output file carries an extension the plotting library cannot
// write.
//
// The suggestion is chosen by edit distance from the supported set, so a
// near miss should be corrected to the extension it most resembles
// rather than to an arbitrary one.
func TestDrawPVSuggestsTheNearestExtension(t *testing.T) {
	ethane, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	cfg := &state.PVConfig{Type: &cubic.PR{}}
	dir := t.TempDir()

	// ".sgv" is a transposition of ".svg" and strictly nearer to it than
	// to any other supported extension. ".pgn" is deliberately omitted:
	// it sits two edits from both ".png" and ".pdf", so which is
	// suggested is a matter of the tie-break rather than of correctness.
	testCases := []struct {
		name    string
		file    string
		suggest string
	}{
		{"transposed svg", "diagram.sgv", ".svg"},
		{"unsupported format", "diagram.txt", ""},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			err := state.DrawPV(cfg, filepath.Join(dir, tc.file), ethane)
			if err == nil {
				t.Fatal("expected an error for an unsupported extension; got nil")
			}

			if !strings.Contains(err.Error(), "Did you mean") {
				t.Fatalf("error = %q; want a suggestion", err)
			}

			if tc.suggest != "" && !strings.Contains(err.Error(), tc.suggest) {
				t.Errorf("error = %q; want it to suggest %q", err, tc.suggest)
			}
		})
	}
}

// TestDrawPVWritesAFile checks that a diagram is produced and is
// non-trivial, for each supported format.
//
// The content is not inspected beyond its opening bytes: the point is
// that the plotting pipeline runs end to end and commits something to
// disk, not that any particular curve was drawn.
func TestDrawPVWritesAFile(t *testing.T) {
	first, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	second, err := state.NewState(substance.Ethane, 490, 70)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	testCases := []struct {
		name   string
		file   string
		prefix string
	}{
		{"svg", "diagram.svg", "<?xml"},
		{"png", "diagram.png", "\x89PNG"},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			output := filepath.Join(t.TempDir(), tc.file)

			cfg := &state.PVConfig{
				Type:           &cubic.PR{},
				Title:          "PV Diagram for Ethane",
				NumberStates:   true,
				LabelIsotherms: true,
			}

			if err := state.DrawPV(cfg, output, first, second); err != nil {
				t.Fatalf("DrawPV returned an unexpected error: %v", err)
			}

			content, err := os.ReadFile(output)
			if err != nil {
				t.Fatalf("the diagram was not written: %v", err)
			}

			if len(content) < 1024 {
				t.Errorf("the diagram is %d bytes, which is too small to contain a plot", len(content))
			}

			if !strings.HasPrefix(string(content), tc.prefix) {
				t.Errorf("the file does not begin as a %s", tc.name)
			}
		})
	}
}

// TestDrawPVDefaultsTheTitle checks that a diagram without a configured
// title is named after its substance.
func TestDrawPVDefaultsTheTitle(t *testing.T) {
	ethane, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	output := filepath.Join(t.TempDir(), "untitled.svg")

	if err := state.DrawPV(&state.PVConfig{Type: &cubic.PR{}}, output, ethane); err != nil {
		t.Fatalf("DrawPV returned an unexpected error: %v", err)
	}

	content, err := os.ReadFile(output)
	if err != nil {
		t.Fatalf("the diagram was not written: %v", err)
	}

	if !strings.Contains(string(content), "Ethane") {
		t.Error("the default title does not name the substance")
	}
}

// TestDrawPVSuggestionIsDeterministic checks that the same typo always
// draws the same suggestion.
//
// Several extensions can be equally near a mistyped one — ".pgn" is two
// edits from both ".png" and ".pdf" — and the candidates are held in a
// map. Ranging over it directly would make the advice depend on the
// iteration order, so the same command could be answered differently on
// consecutive runs.
func TestDrawPVSuggestionIsDeterministic(t *testing.T) {
	ethane, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		t.Fatalf("NewState returned an unexpected error: %v", err)
	}

	cfg := &state.PVConfig{Type: &cubic.PR{}}
	output := filepath.Join(t.TempDir(), "diagram.pgn")

	first := state.DrawPV(cfg, output, ethane)
	if first == nil {
		t.Fatal("expected an error for an unsupported extension; got nil")
	}

	for range 20 {
		again := state.DrawPV(cfg, output, ethane)
		if again == nil {
			t.Fatal("expected an error; got nil")
		}

		if again.Error() != first.Error() {
			t.Fatalf("the suggestion varies between runs: %v then %v", first, again)
		}
	}
}
