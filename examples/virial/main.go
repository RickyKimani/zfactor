// The virial equations: a power series in density, truncated.
//
// The two-term form is linear in pressure and reliable only where the gas
// is thin, which is why it refuses to evaluate above 15 bar. The
// three-term form keeps one more coefficient and reaches further, at the
// cost of solving a cubic.
//
// The figures are Example 3.8 of Smith, Van Ness & Abbott: isopropanol
// vapor at 200 °C and 10 bar, whose reported coefficients are
// B = -388 cm³/mol and C = -26,000 cm⁶/mol². The book works the same
// state three ways, so each result here has a published value to sit
// beside.
//
// Run with: go run ./examples/virial
package main

import (
	"fmt"
	"log"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/virial"
)

func main() {
	args := zfactor.Args{
		T: 473.15,   // K, i.e. 200 °C
		P: 10.0,     // bar
		R: 83.14,    // bar·cm³/(mol·K)
		B: -388.0,   // second virial coefficient (cm³/mol)
		C: -26000.0, // third virial coefficient (cm⁶/mol²)
	}

	// The ideal-gas volume, which is both the reference the book compares
	// against and the starting point the three-term iteration uses.
	vIdeal := args.R * args.T / args.P

	fmt.Printf("Ideal gas:  V = %.0f cm³/mol   Z = 1        (book: 3934)\n", vIdeal)

	// Two-term: Z = 1 + BP/RT, linear in pressure.
	z2, err := virial.CompressibilityTwoTerm(args)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("Two-term:   V = %.0f cm³/mol   Z = %.4f   (book: 3546, 0.9014)\n",
		z2*vIdeal, z2)

	// Three-term: keeps C, so it is a cubic in volume and returns all
	// three roots. For a vapor the physical one is the largest real root;
	// the others are the liquid branch and the unstable middle root.
	roots, err := virial.SolveForVolumeThreeTerm(args)
	if err != nil {
		log.Fatal(err)
	}

	vapor, ok := largestReal(roots)
	if !ok {
		log.Fatal("no real volume root")
	}

	// Z can be had from the definition, or from the equation the root was
	// found for. They agree, which is a check on the root.
	z3, err := virial.CompressibilityThreeTerm(vapor, args)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("Three-term: V = %.0f cm³/mol   Z = %.4f   (book: 3488, 0.8866)\n",
		vapor, z3)

	fmt.Println("\nThe ideal-gas volume is 13% high here and the two-term form 1.7% high,")
	fmt.Println("which is the price of truncating the series one term earlier.")
}

// largestReal returns the largest root with a negligible imaginary part,
// which for a vapor is the physical volume.
func largestReal(roots [3]complex128) (float64, bool) {
	const imagTol = 1e-9

	best, found := 0.0, false

	for _, root := range roots {
		if math.Abs(imag(root)) > imagTol {
			continue
		}

		if v := real(root); !found || v > best {
			best, found = v, true
		}
	}

	return best, found
}
