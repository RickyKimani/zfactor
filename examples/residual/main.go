// Residual properties: how far a real fluid departs from an ideal gas.
//
// Three routes to the same two quantities are compared. The Abbott
// correlations derive them from the second virial coefficient and are adequate
// at low reduced pressure; the Lee-Kesler tables are needed above it; and a
// cubic equation of state computes them from its own parameters, carrying no
// tables at all.
//
// All three return the same dimensionless groups, H^R/(R Tc) and S^R/R, so
// they can be set side by side.
//
// Run with: go run ./examples/residual
package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
	"github.com/rickykimani/zfactor/substance"
)

// R in bar*cm^3/(mol*K), which the equations of state need because they solve
// for a volume. The correlations do not: they work in reduced coordinates.
const R = 10 * zfactor.RSI

func main() {
	correlations()
	fmt.Println()

	fromEquationsOfState()
	fmt.Println()

	phaseMatters()
}

// correlations shows the two generalized routes, which need only the reduced
// temperature and pressure.
func correlations() {
	eth := substance.Ethane
	args := zfactor.Args{T: 299.0, P: 32.0}

	fmt.Println("Ethane at 299 K and 32 bar, from generalized correlations")

	hR, err := eth.AbbottResidualEnthalpy(args)
	if err != nil {
		log.Fatal(err)
	}

	sR, err := eth.AbbottResidualEntropy(args)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("  Abbott       H^R/(R Tc) = %8.4f    S^R/R = %8.4f\n", hR, sR)

	hrLK, err := eth.LeeKesler(args, leekesler.ResidualEnthalpy)
	if err != nil {
		log.Fatal(err)
	}

	srLK, err := eth.LeeKesler(args, leekesler.ResidualEntropy)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("  Lee-Kesler   H^R/(R Tc) = %8.4f    S^R/R = %8.4f\n", hrLK, srLK)
	fmt.Println("\n  Abbott is a second-virial correlation, so it is the one to")
	fmt.Println("  distrust first as the pressure rises.")
}

// fromEquationsOfState reproduces Table 13.7 of Smith, Van Ness and Abbott:
// n-butane at 500 K and 50 bar, worked for each of the four equations.
//
// The published values are printed alongside, which is the useful form for an
// example: it shows both how to call the functions and how far the equations
// disagree with each other and with measurement.
func fromEquationsOfState() {
	butane := substance.NButane

	// The book's gas constant, in J/(mol*K), for converting the
	// dimensionless groups into an enthalpy and an entropy.
	const bookR = 8.314

	Tc := butane.Critical.Tc
	args := zfactor.Args{T: 500.0, P: 50.0, R: R}

	equations := []struct {
		name string
		eos  cubic.EOSType
		// From Table 13.7, in J/mol and J/(mol*K).
		bookH float64
		bookS float64
	}{
		{"vdW", &cubic.VdW{}, -3937, -5.424},
		{"RK", &cubic.RK{}, -4505, -6.546},
		{"SRK", &cubic.SRK{}, -4824, -7.413},
		{"PR", &cubic.PR{}, -4988, -7.426},
	}

	fmt.Println("n-Butane at 500 K and 50 bar, from equations of state")
	fmt.Println("  (Table 13.7 of Smith, Van Ness & Abbott)")
	fmt.Println()
	fmt.Printf("  %-6s %12s %12s   %12s %12s\n", "", "H^R J/mol", "book", "S^R J/mol/K", "book")

	for _, e := range equations {
		hr, err := butane.CubicResidualEnthalpy(e.eos, cubic.StablePhase, args)
		if err != nil {
			log.Fatal(err)
		}

		sr, err := butane.CubicResidualEntropy(e.eos, cubic.StablePhase, args)
		if err != nil {
			log.Fatal(err)
		}

		// H^R/(R Tc) back to an enthalpy, and S^R/R to an entropy.
		fmt.Printf("  %-6s %12.1f %12.1f   %12.3f %12.3f\n",
			e.name, hr*bookR*Tc, e.bookH, sr*bookR, e.bookS)
	}

	// The residual Gibbs energy is the third departure property, and it ties
	// the other two together: ln phi = H^R/RT - S^R/R. Checking that identity
	// is a cheap way to confirm a state was set up as intended.
	srk := &cubic.SRK{}

	hr, err := butane.CubicResidualEnthalpy(srk, cubic.StablePhase, args)
	if err != nil {
		log.Fatal(err)
	}

	sr, err := butane.CubicResidualEntropy(srk, cubic.StablePhase, args)
	if err != nil {
		log.Fatal(err)
	}

	lnPhi, err := butane.CubicLogFugacity(srk, cubic.StablePhase, args)
	if err != nil {
		log.Fatal(err)
	}

	// The enthalpy is normalised by R Tc and the identity is in terms of
	// H^R/(R T), so the reduced temperature comes back out.
	overRT := hr / (args.T / Tc)

	fmt.Printf("\n  SRK: ln phi = %.6f, and H^R/RT - S^R/R = %.6f\n", lnPhi, overRT-sr)
}

// phaseMatters shows why the phase is an argument rather than a default.
//
// Below the critical temperature the equation has a liquid root and a vapour
// root, and each has its own residual properties. At the saturation pressure
// both exist at once, which is the case a single answer cannot describe.
func phaseMatters() {
	butane := substance.NButane

	T := 0.85 * butane.Critical.Tc

	// SaturationPressure is given the temperature and finds the pressure, so
	// the configuration it takes carries no pressure of its own.
	eos := &cubic.PR{}
	cfg := butane.CubicConfig(eos, zfactor.Args{T: T, R: R})

	pSat, err := cubic.SaturationPressure(cfg, T)
	if err != nil {
		log.Fatal(err)
	}

	args := zfactor.Args{T: T, P: pSat, R: R}

	fmt.Printf("n-Butane at %.1f K, on the saturation line at %.2f bar (Peng-Robinson)\n", T, pSat)

	for _, phase := range []cubic.Phase{cubic.LiquidPhase, cubic.VaporPhase} {
		hr, err := butane.CubicResidualEnthalpy(eos, phase, args)
		if err != nil {
			log.Fatal(err)
		}

		sr, err := butane.CubicResidualEntropy(eos, phase, args)
		if err != nil {
			log.Fatal(err)
		}

		fmt.Printf("  %-6s  H^R/(R Tc) = %8.4f    S^R/R = %8.4f\n", phase, hr, sr)
	}

	fmt.Println("\n  The gap between the two is the enthalpy of vaporisation, which is")
	fmt.Println("  why a single answer for this state would have to be wrong for one")
	fmt.Println("  phase. cubic.StablePhase picks whichever is stable by comparing")
	fmt.Println("  fugacities, and on the saturation line those are equal.")
}
