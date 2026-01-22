package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/substance"
)

func main() {
	// Saturation Pressure (Antoine Equation)
	// Note: Antoine coefficients often use Celsius and specific pressure units (e.g., kPa)
	// Here we calculate saturation pressure of Ethanol at 25°C
	pSat, err := antoine.Ethanol.Pressure(25.0) // 25°C
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Saturation Pressure (Ethanol @ 25C): %.2f kPa\n", pSat)

	// Saturated Liquid Volume (Rackett Equation)
	eth := substance.Ethane
	vSat, err := eth.Vsat(299.0) // T in Kelvin required
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Saturated Liquid Volume (Ethane @ 299K): %.4f cm³/mol\n", vSat)

	// Reduced Density (Lydersen Charts)
	rhoR, err := eth.ReducedDensity(zfactor.Args{T: 299.0, P: 50.0})
	if err != nil {
		log.Fatal(err)
	}
	fmt.Printf("Reduced Density (Ethane @ 299K, 50bar): %.4f\n", rhoR)
}
