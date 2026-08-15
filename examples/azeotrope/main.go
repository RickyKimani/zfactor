// Finding azeotropes: the compositions at which a mixture boils without
// changing composition.
//
// At an azeotrope the liquid and vapor compositions coincide, so
// distillation cannot separate the mixture past that point. Locating them
// is what tells you where a separation will stall.
//
// Run with: go run ./examples/azeotrope
package main

import (
	"errors"
	"fmt"
	"log"

	"github.com/rickykimani/zfactor/antoine"
	modified_raoult "github.com/rickykimani/zfactor/vle/modified-raoult"
)

// The methanol(1)/methyl acetate(2) system of Example 13.1, whose Antoine
// constants the book gives with temperature in kelvin. This library works
// in °C, and since the temperature appears only as (T + C) the conversion
// is absorbed entirely into C: A and B are unchanged.
func exampleModels() []antoine.Model {
	return []antoine.Model{
		&antoine.Antoine{
			Name: "methanol",
			A:    16.59158, B: 3643.31, C: 273.15 - 33.424,
			Range: antoine.TempRange{Low: 0, High: 150},
		},
		&antoine.Antoine{
			Name: "methyl acetate",
			A:    14.25326, B: 2665.54, C: 273.15 - 53.424,
			Range: antoine.TempRange{Low: 0, High: 150},
		},
	}
}

func main() {
	isothermal()
	fmt.Println()

	isobaric()
	fmt.Println()

	noAzeotrope()
	fmt.Println()

	twoAzeotropes()
}

// isothermal finds the azeotrope at a fixed temperature, part (e) of
// Example 13.1.
func isothermal() {
	// The example's Margules parameter varies with temperature; at
	// 318.15 K, which is 45 °C, it is:
	const a = 2.771 - 0.00523*318.15

	found, err := modified_raoult.AzeotropeP(modified_raoult.MixtureInput{
		T:        45,
		Antoine:  exampleModels(),
		Activity: modified_raoult.Margules{A12: a, A21: a},
	})
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("At a fixed temperature — methanol/methyl acetate at 45 °C")

	for _, azeotrope := range found {
		fmt.Printf(
			"  x1 = y1 = %.4f at %.2f kPa\n",
			azeotrope.X[0], azeotrope.P,
		)
	}

	fmt.Println("\n  the book gives x1 = 0.325 at 73.8 kPa")
}

// isobaric finds the azeotrope at a fixed pressure, where the temperature
// is what is being solved for.
//
// Each trial composition needs a bubble-temperature calculation before
// the azeotropic condition can be evaluated, so this costs more than the
// isothermal case.
func isobaric() {
	const a = 2.771 - 0.00523*318.15

	found, err := modified_raoult.AzeotropeT(modified_raoult.MixtureInput{
		P:        73.8,
		Antoine:  exampleModels(),
		Activity: modified_raoult.Margules{A12: a, A21: a},
	})
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("At a fixed pressure — the same system at 73.8 kPa")

	for _, azeotrope := range found {
		fmt.Printf(
			"  x1 = y1 = %.4f at %.2f °C\n",
			azeotrope.X[0], azeotrope.T,
		)
	}
}

// noAzeotrope shows the common outcome: most mixtures do not form one.
//
// It is reported as a sentinel error rather than an empty result, so a
// caller can distinguish "no azeotrope exists" from a failure to converge.
func noAzeotrope() {
	// An ideal solution has unit activity coefficients, so its relative
	// volatility never reaches one and no azeotrope can form.
	_, err := modified_raoult.AzeotropeP(modified_raoult.MixtureInput{
		T:        45,
		Antoine:  exampleModels(),
		Activity: modified_raoult.Margules{A12: 0, A21: 0},
	})

	fmt.Println("An ideal solution")

	if errors.Is(err, modified_raoult.ErrNoAzeotrope) {
		fmt.Printf("  %v\n", err)
		fmt.Println("  — as expected: without activity coefficients the volatility")
		fmt.Println("    ratio is fixed by the saturation pressures alone.")
		return
	}

	log.Fatalf("expected no azeotrope, got %v", err)
}

// twoAzeotropes shows why the result is a slice.
//
// A binary can form more than one, and a search that only compared the
// two pure limits would miss both: with an even number of azeotropes the
// residual has the same sign at each end. The parameters here are chosen
// to produce a pair.
func twoAzeotropes() {
	found, err := modified_raoult.AzeotropeP(modified_raoult.SaturationPressureInput{
		T:        45,
		PSats:    []float64{30.0, 88.35},
		Activity: modified_raoult.Margules{A12: 1.0, A21: 3.0},
	})
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("A double azeotrope — asymmetric Margules parameters")
	fmt.Printf("  %d azeotropes found\n", len(found))

	for i, azeotrope := range found {
		fmt.Printf(
			"    %d: x1 = y1 = %.6f at %.4f kPa\n",
			i+1, azeotrope.X[0], azeotrope.P,
		)
	}

	fmt.Println("\n  Both lie strictly inside (0, 1), and the residual has the same")
	fmt.Println("  sign at each pure limit — which is why the search sweeps the")
	fmt.Println("  composition range rather than testing the endpoints.")
}
