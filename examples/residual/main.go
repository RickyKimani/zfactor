package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor"
	leekesler "github.com/rickykimani/zfactor/lee-kesler"
	"github.com/rickykimani/zfactor/substance"
)

func main() {
	eth := substance.Ethane
	args := zfactor.Args{T: 299.0, P: 32.0}

	fmt.Println("--- Abbott (Virial) correlations ---")
	// Dimensionless Residual Enthalpy (H^R / RTc)
	// Uses Abbott/Virial generalized correlations
	hR, err := eth.ResidualEnthalpy(args)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Residual Enthalpy (H^R / RTc): %.4f\n", hR)

	// Dimensionless Residual Entropy (S^R / R)
	sR, err := eth.ResidualEntropy(args)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Residual Entropy (S^R / R): %.4f\n", sR)

	fmt.Println("\n--- Lee-Kesler correlations ---")
	// Lee-Kesler generally provides more accurate values at higher pressures

	// Dimensionless Residual Enthalpy (H^R / RTc)
	hR_LK, err := eth.LeeKesler(args, leekesler.ResidualEnthalpy)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Residual Enthalpy (H^R / RTc): %.4f\n", hR_LK)

	// Dimensionless Residual Entropy (S^R / R)
	sR_LK, err := eth.LeeKesler(args, leekesler.ResidualEntropy)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Residual Entropy (S^R / R): %.4f\n", sR_LK)

}
