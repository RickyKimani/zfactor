package main

import (
	"fmt"
	"log"

	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/vle/raoult"
)

func main() {

	mi := raoult.MixtureInput{
		T:            100,
		P:            120,
		Compositions: []float64{0.33, 1 - 0.33},
		Antoine:      []antoine.Model{antoine.Benzene, antoine.Toluene},
	}
	bpr, err := raoult.BubbleP(mi)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(bpr)

	btr, err := raoult.BubbleT(mi)
	if err != nil {
		log.Fatal(err)
	}
	fmt.Println(btr)

}
