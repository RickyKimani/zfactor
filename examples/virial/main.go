package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/virial"
)

func main() {
	// Isopropanol vapor example
	args := zfactor.Args{
		T: 473.15,   // K
		P: 10.0,     // bar
		R: 83.14,    // bar·cm³/(mol·K)
		B: -338.0,   // Second virial coefficient (cm³/mol)
		C: -26000.0, // Third virial coefficient (cm⁶/mol²)
	}

	// 2-Term Virial (Z = 1 + BP/RT)
	// Used for low to moderate pressures
	z2, err := virial.CompressibilityTwoTerm(args)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Z (2-term): %.4f\n", z2)

	// 3-Term Virial (Iterative solution)
	// Returns complex roots for volume
	// Used for higher pressures
	roots, err := virial.SolveForVolumeThreeTerm(args)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Volume Roots (3-term): %v\n", roots)
}
