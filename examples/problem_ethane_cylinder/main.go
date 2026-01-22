package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
	"github.com/rickykimani/zfactor/state"
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
		P1 = 32               // bar
		T1 = 299              // K
		R  = 10 * zfactor.RSI // bar*cm³/(mol*K)
		T2 = 490              // K
	)

	// Initialize the initial thermodynamic state (State 1) at 299 K and 32 bar.
	s1, err := state.NewState(ethane, T1, P1)
	if err != nil {
		log.Fatal(err)
	}

	// Compute the compressibility factor (Z) using the Lee-Kesler correlation.
	// This method is suitable here as the state is in the single-phase region.
	z, err := ethane.LeeKesler(zfactor.Args{T: s1.Temperature, P: s1.Pressure}, leekesler.CompressibilityFactor)
	if err != nil {
		log.Fatal(err)
	}

	// Calculate the molar volume (v) using the definition of Z (v = ZRT/P).
	// Since the system is a closed cylinder, the process is isochoric (constant volume), so v1 = v2.
	v := z * R * T1 / P1
	fmt.Printf("Molar Volume at T1: %.4f cm³/mol\n", v)

	// Configure the Soave-Redlich-Kwong (SRK) Equation of State for the final temperature (T2).
	// Pressure is initialized to 0 as it is the variable to be determined.
	cfg := ethane.CubicConfig(&cubic.SRK{}, zfactor.Args{T: T2, R: R})

	// Calculate the final pressure (P2) corresponding to the constant molar volume using the SRK EOS.
	// The result includes the calculated pressure and intermediate EOS parameters.
	pressureResult, err := cubic.Pressure(cfg, v)
	if err != nil {
		log.Fatal(err)
	}

	// Retrieve the calculated final pressure.
	P2 := pressureResult.P

	// Output the calculated final pressure.
	fmt.Printf("Final Pressure P2 = %.4f bar\n", P2)

	// Initialize the final thermodynamic state (State 2) using the final temperature and calculated pressure.
	s2, err := state.NewState(ethane, T2, P2)
	if err != nil {
		log.Fatal(err)
	}

	// Configure the visualization parameters for the PV diagram, enabling state numbering.
	pvCfg := &state.PVConfig{
		Type:                  cfg.Type,
		NumberStates:          true,
		LabelIsotherms:        true,
		TitleColor:            state.Green,
		IsothermsColor:        state.Purple,
		IsothermLabelColor:    state.Orange,
		CriticalIsothermColor: state.Red,
		StatePointColor:       state.Blue,
		ShowOutputPath:        true,
	}

	// Generate and save the PV diagram to the specified output file.
	err = state.DrawPV(pvCfg, "pv.png", s1, s2)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println("Generated PV diagram at pv.png")
}
