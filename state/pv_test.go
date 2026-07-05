package state_test

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/state"
	"github.com/rickykimani/zfactor/state/palettes"
	"github.com/rickykimani/zfactor/state/themes"
	"github.com/rickykimani/zfactor/substance"
)

func TestDrawPVThemes(t *testing.T) {
	sub := substance.Ethane

	// Create some states for Ethane (Tc = 305.3 K, Pc = 48.72 bar)
	s1, err := state.NewState(sub, 280.0, 20.0) // subcritical state
	if err != nil {
		t.Fatalf("failed to create state 1: %v", err)
	}
	s2, err := state.NewState(sub, 320.0, 60.0) // supercritical state
	if err != nil {
		t.Fatalf("failed to create state 2: %v", err)
	}

	eos := &cubic.PR{}

	// Test default theme (light theme)
	cfgLight := state.DefaultPVConfig(eos)
	cfgLight.Theme = themes.Default
	cfgLight.NumberStates = true
	cfgLight.LabelIsotherms = true
	cfgLight.ShowOutputPath = false

	outputPathLight := filepath.Join(t.TempDir(), "pv_light.svg")

	err = state.DrawPV(cfgLight, outputPathLight, s1, s2)
	if err != nil {
		t.Fatalf("DrawPV failed with Default theme: %v", err)
	}

	if _, err := os.Stat(outputPathLight); os.IsNotExist(err) {
		t.Errorf("expected file %s to be created, but it does not exist", outputPathLight)
	}

	// Test dark theme with Viridis palette
	cfgDark := state.DefaultPVConfig(eos)
	cfgDark.Theme = themes.Dark.WithPalette(palettes.Viridis)
	cfgDark.NumberStates = true
	cfgDark.LabelIsotherms = true
	cfgDark.ShowOutputPath = false

	outputPathDark := filepath.Join(t.TempDir(), "pv_dark.svg")

	err = state.DrawPV(cfgDark, outputPathDark, s1, s2)
	if err != nil {
		t.Fatalf("DrawPV failed with Dark theme: %v", err)
	}

	if _, err := os.Stat(outputPathDark); os.IsNotExist(err) {
		t.Errorf("expected file %s to be created, but it does not exist", outputPathDark)
	}
}
