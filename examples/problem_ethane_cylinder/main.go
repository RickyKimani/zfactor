package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
	"github.com/rickykimani/zfactor/state"
	"github.com/rickykimani/zfactor/state/palettes"
	"github.com/rickykimani/zfactor/state/themes"
	"github.com/rickykimani/zfactor/substance"
)

/*
Problem Statement:
The vapor pressure of ethane at 299 K is 42.7 bar.
A closed cylinder contains ethane at 299 K and 32 bar.
The cylinder is subsequently heated to 490 K.

Tasks:
a) Identify the thermodynamic states of ethane at 299 K and 490 K using a PV diagram.
b) Determine the molar volume of ethane in the cylinder at 299 K.
c) Determine the pressure of ethane in the cylinder at 490 K.
*/

func main() {
	ethane := substance.Ethane

	const (
		P1 = 32.0             // bar
		T1 = 299.0            // K
		T2 = 490.0            // K
		R  = 10 * zfactor.RSI // bar·cm³/(mol·K)
	)

	// ------------------------------------------------------------
	// Initial state
	// ------------------------------------------------------------

	s1, err := state.NewState(ethane, T1, P1)
	if err != nil {
		log.Fatal(err)
	}

	// Compute the compressibility factor (Z) using the Lee-Kesler correlation.
	// This method is suitable here as the state is in the single-phase region.
	z, err := ethane.LeeKesler(
		zfactor.Args{
			T: s1.Temperature,
			P: s1.Pressure,
		},
		leekesler.CompressibilityFactor,
	)
	if err != nil {
		log.Fatal(err)
	}

	// Calculate the molar volume (v) using the definition of Z (v = ZRT/P).
	// Since the system is a closed cylinder, the process is isochoric (constant volume), so v1 = v2.
	v := z * R * T1 / P1
	fmt.Printf("Molar volume at %.0f K = %.4f cm³/mol\n", T1, v)

	// ------------------------------------------------------------
	// Final state
	// ------------------------------------------------------------

	cfg := ethane.CubicConfig(
		&cubic.SRK{},
		zfactor.Args{
			T: T2,
			R: R,
		},
	)

	pressureResult, err := cubic.Pressure(cfg, v)
	if err != nil {
		log.Fatal(err)
	}

	P2 := pressureResult.P
	fmt.Printf("Pressure at %.0f K = %.4f bar\n", T2, P2)

	s2, err := state.NewState(ethane, T2, P2)
	if err != nil {
		log.Fatal(err)
	}

	// ------------------------------------------------------------
	// Draw PV diagram
	// ------------------------------------------------------------

	pvCfg := state.DefaultPVConfig(cfg.Type)

	pvCfg.NumberStates = true
	pvCfg.LabelIsotherms = true

	// Optional styling
	pvCfg.Theme = themes.
		Default.
		WithPalette(palettes.Viridis)

	err = state.DrawPV(pvCfg, "pv.png", s1, s2)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Generated PV diagram at pv.png")
}
