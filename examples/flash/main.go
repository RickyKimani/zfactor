// Flash calculations: what a feed separates into at a given temperature
// and pressure.
//
// Bubble- and dew-point calculations locate the edges of the two-phase
// region. A flash answers the question between them — how much of the
// feed is vapor, and what each phase contains.
//
// Run with: go run ./examples/flash
package main

import (
	"errors"
	"fmt"
	"log"

	"github.com/rickykimani/zfactor/antoine"
	"github.com/rickykimani/zfactor/vle"
	modified_raoult "github.com/rickykimani/zfactor/vle/modified-raoult"
	"github.com/rickykimani/zfactor/vle/raoult"
)

func main() {
	idealFlash()
	fmt.Println()

	theTwoPhaseWindow()
	fmt.Println()

	nonIdealFlash()
}

// idealFlash reproduces Example 13.8 of Smith, Van Ness & Abbott:
// acetone/acetonitrile/nitromethane at 80 °C and 110 kPa.
//
// The saturation pressures are supplied directly rather than through a
// correlation, which is what SaturationPressureInput is for.
func idealFlash() {
	fmt.Println("Raoult's law flash — acetone(1)/acetonitrile(2)/nitromethane(3)")
	fmt.Println("  80 °C, 110 kPa, feed 0.45 / 0.35 / 0.20")

	result, err := raoult.FlashPT(raoult.SaturationPressureInput{
		P:            110,
		Compositions: []float64{0.45, 0.35, 0.20},
		PSats:        []float64{195.75, 97.84, 50.32},
	})
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("\n  vapor  %.4f mol   liquid %.4f mol\n", result.V, result.L)
	fmt.Println("\n  species          x        y")

	names := []string{"acetone", "acetonitrile", "nitromethane"}
	for i, name := range names {
		fmt.Printf("  %-14s %.4f   %.4f\n", name, result.X[i], result.Y[i])
	}

	fmt.Println("\n  the book gives V = 0.7364, x = 0.2859/0.3810/0.3331,")
	fmt.Println("                          y = 0.5087/0.3389/0.1524")
}

// theTwoPhaseWindow shows what happens either side of the region where a
// feed separates.
//
// A flash has no answer outside it, and the error says which side you are
// on so that a caller adjusting conditions knows which way to move.
func theTwoPhaseWindow() {
	feed := []float64{0.45, 0.35, 0.20}
	psat := []float64{195.75, 97.84, 50.32}

	// The boundaries: the feed as a liquid gives the bubble pressure, the
	// same composition as a vapor gives the dew pressure.
	bubble, err := raoult.BubbleP(raoult.SaturationPressureInput{
		Compositions: feed, PSats: psat,
	})
	if err != nil {
		log.Fatal(err)
	}

	dew, err := raoult.DewP(raoult.SaturationPressureInput{
		Compositions: feed, PSats: psat,
	})
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("The feed separates between %.2f and %.2f kPa\n", dew.P, bubble.P)
	fmt.Println("\n  pressure   vapor fraction")

	for _, pressure := range []float64{60, 95, 105, 110, 125, 140} {
		result, err := raoult.FlashPT(raoult.SaturationPressureInput{
			P: pressure, Compositions: feed, PSats: psat,
		})

		switch {
		case err == nil:
			fmt.Printf("  %6.1f     %.4f\n", pressure, result.V)

		default:
			// A feed that does not split is described rather than
			// reported as a failure.
			var single *vle.SinglePhaseError
			if errors.As(err, &single) {
				fmt.Printf("  %6.1f     — %s\n", pressure, single.State)
				continue
			}

			log.Fatal(err)
		}
	}
}

// nonIdealFlash flashes a mixture whose liquid phase is not ideal.
//
// The equilibrium ratios now depend on the liquid composition being
// solved for, so the calculation iterates: a bubble- and a dew-point
// calculation establish the range and the starting point, then the vapor
// fraction and the activity coefficients are refined together.
func nonIdealFlash() {
	// Methanol/methyl acetate at 45 °C, the system of Example 13.1. Its
	// Margules parameter varies with temperature; at 318.15 K it is:
	const a = 2.771 - 0.00523*318.15

	feed := []float64{0.40, 0.60}

	input := modified_raoult.MixtureInput{
		T:            45,
		Compositions: feed,
		Antoine: []antoine.Model{
			antoine.Methanol,
			antoine.MethylAcetate,
		},
		Activity: modified_raoult.Margules{A12: a, A21: a},
	}

	bubble, err := modified_raoult.BubbleP(input)
	if err != nil {
		log.Fatal(err)
	}

	dew, err := modified_raoult.DewP(input)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Modified Raoult's law flash — methanol(1)/methyl acetate(2)")
	fmt.Printf("  45 °C, feed %.2f / %.2f, Margules A = %.5f\n", feed[0], feed[1], a)
	fmt.Printf("  separates between %.2f and %.2f kPa\n", dew.P, bubble.P)
	fmt.Println("\n  pressure   vapor      x1       y1")

	// Sample across the window.
	for _, fraction := range []float64{0.1, 0.3, 0.5, 0.7, 0.9} {
		pressure := dew.P + fraction*(bubble.P-dew.P)

		flashed := input
		flashed.P = pressure

		result, err := modified_raoult.FlashPT(flashed)
		if err != nil {
			log.Fatal(err)
		}

		fmt.Printf(
			"  %7.2f    %.4f   %.4f   %.4f\n",
			pressure, result.V, result.X[0], result.Y[0],
		)
	}

	fmt.Println("\n  Note the limits: at the dew pressure almost everything is")
	fmt.Println("  vapor and the liquid is the heavier component; at the bubble")
	fmt.Println("  pressure the liquid is the feed.")
}
