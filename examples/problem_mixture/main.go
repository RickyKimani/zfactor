package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
	"github.com/rickykimani/zfactor/substance"
)

func main() {
	// Estimate V, H^R and S^R for an equimolar mixture
	// of carbon dioxide and propane at 450K and 140 bar by the Lee/Kesler correlations.

	args := zfactor.Args{T: 450, P: 140}
	components := []substance.Component{
		{Substance: substance.CarbonDioxide, Fraction: 0.5},
		{Substance: substance.Propane, Fraction: 0.5},
	}
	mixture, err := substance.NewLinearMixture("mixture", components)
	if err != nil {
		log.Fatal(err)
	}

	// Calculate properties
	fmt.Printf("Pseudo-Critical Properties: Tc=%.2f K, Pc=%.2f bar, Omega=%.4f\n",
		mixture.Critical.Tc, mixture.Critical.Pc, mixture.Acentric)
	fmt.Println("------------------------------------------------")

	// 1. Compressibility Factor (Z)
	// We use the pseudo-critical properties to find Z from Lee-Kesler tables.
	z, err := mixture.LeeKesler(args, leekesler.CompressibilityFactor)
	if err != nil {
		log.Fatal(err)
	}

	// Calculate Molar Volume: V = ZRT/P
	// R = 83.14 bar*cm³/(mol*K) (using zfactor.RSI * 10)
	R := zfactor.RSI * 10
	v := z * R * args.T / args.P
	fmt.Printf("Compressibility Factor (Z): %.4f\n", z)
	fmt.Printf("Molar Volume (V):           %.4f cm³/mol\n", v)

	// 2. Residual Enthalpy (H^R)
	// Lee-Kesler returns dimensionless H^R / (R * Tc_mix)
	h_dim, err := mixture.LeeKesler(args, leekesler.ResidualEnthalpy)
	if err != nil {
		log.Fatal(err)
	}

	// Convert to dimensional form: H^R = (H^R/RTc) * R * Tc
	HR := h_dim * zfactor.RSI * mixture.Critical.Tc
	fmt.Printf("Residual Enthalpy (H^R):    %.2f J/mol\n", HR)

	// 3. Residual Entropy (S^R)
	// Lee-Kesler returns dimensionless S^R / R
	s_dim, err := mixture.LeeKesler(args, leekesler.ResidualEntropy)
	if err != nil {
		log.Fatal(err)
	}

	// Convert to dimensional form: S^R = (S^R/R) * R
	SR := s_dim * zfactor.RSI
	fmt.Printf("Residual Entropy (S^R):     %.2f J/(mol·K)\n", SR)
}
