package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor/antoine"
	modified_raoult "github.com/rickykimani/zfactor/vle/modified-raoult"
	"github.com/rickykimani/zfactor/vle/raoult"
)

func main() {
	// ------------------------------------------------------------
	// Raoult's law example
	// ------------------------------------------------------------

	x := []float64{0.30, 0.70}

	ideal := raoult.MixtureInput{
		T:            100,
		P:            101.33,
		Compositions: x,
		Antoine: []antoine.Model{
			antoine.Benzene,
			antoine.Toluene,
		},
	}

	bp, err := raoult.BubbleP(ideal)
	if err != nil {
		log.Fatal(err)
	}

	bt, err := raoult.BubbleT(ideal)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Raoult Bubble P:", bp)
	fmt.Println("Raoult Bubble T:", bt)

	// ------------------------------------------------------------
	// Modified Raoult's law (Wilson activity model)
	// ------------------------------------------------------------

	molarVolumes := []float64{
		74.05,
		18.07,
	}

	// Wilson interaction parameters (cal/mol)
	aCal := [][]float64{
		{0, 291.27},
		{1448.01, 0},
	}

	const calToJ = 4.186

	a := make([][]float64, len(aCal))
	for i := range aCal {
		a[i] = make([]float64, len(aCal[i]))
		for j := range aCal[i] {
			a[i][j] = aCal[i][j] * calToJ
		}
	}

	nonIdeal := modified_raoult.MixtureInput{
		P:            101.33,
		Compositions: x,
		Antoine: []antoine.Model{
			antoine.Acetone,
			antoine.Water,
		},
		Activity: modified_raoult.Wilson{
			V:           molarVolumes,
			Interaction: a,
		},
	}

	btm, err := modified_raoult.BubbleT(nonIdeal)
	if err != nil {
		log.Fatal(err)
	}

	dtm, err := modified_raoult.DewT(nonIdeal)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Modified Bubble T:", btm)
	fmt.Println("Modified Dew T:", dtm)
}
