// Drawing a PV diagram: the saturation dome, the critical isotherm and
// whatever states you want marked on it.
//
// This is the code the README's example diagram was produced from.
//
// Run with: go run ./examples/pvdiagram
package main

import (
	"fmt"
	"log"
	"path/filepath"

	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/state"
	"github.com/rickykimani/zfactor/state/themes"
	"github.com/rickykimani/zfactor/substance"
)

func main() {
	// Two states of ethane: one inside the two-phase region and one
	// supercritical.
	first, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		log.Fatal(err)
	}

	second, err := state.NewState(substance.Ethane, 490, 70)
	if err != nil {
		log.Fatal(err)
	}

	// Every state on one diagram must be the same substance: the axes are
	// scaled from its critical properties and the dome drawn from its
	// equation of state.
	cfg := &state.PVConfig{
		Type:           &cubic.PR{},
		Title:          "PV Diagram for Ethane",
		NumberStates:   true,
		LabelIsotherms: true,
	}

	// The diagrams are written beside this example. The extension chooses
	// the format; an unsupported one is refused with the nearest
	// supported spelling suggested.
	const dir = "examples/pvdiagram"

	light := filepath.Join(dir, "ethane_pv.svg")

	if err := state.DrawPV(cfg, light, first, second); err != nil {
		log.Fatal(err)
	}

	fmt.Println("wrote", light)

	// The same diagram against a dark background. A theme governs the
	// background, text and axes; the colours of the thermodynamic series
	// come from a palette instead. The two are deliberately independent,
	// so either can be changed without the other.
	dark := *cfg
	dark.Theme = themes.DarkTheme()

	darkFile := filepath.Join(dir, "ethane_pv_dark.svg")

	if err := state.DrawPV(&dark, darkFile, first, second); err != nil {
		log.Fatal(err)
	}

	fmt.Println("wrote", darkFile)

	// A raster format, for somewhere that cannot render SVG.
	raster := filepath.Join(dir, "ethane_pv.png")

	if err := state.DrawPV(cfg, raster, first, second); err != nil {
		log.Fatal(err)
	}

	fmt.Println("wrote", raster)

	fmt.Println("\nThe dome is drawn from the Peng/Robinson equation, so it is a")
	fmt.Println("prediction rather than measured data. Swapping the Type field for")
	fmt.Println("&cubic.VdW{} shows how much the choice of equation moves it.")
}
