// Solving linear systems, and the iteration pattern the rest of the
// library is built to support.
//
// The package exists for the inner loop of an equilibrium or reaction
// calculation, where the same small system is solved over and over. The
// second half of this example is that loop, and the reason the API has
// both Factorize and Refactorize.
//
// Run with: go run ./examples/linalglite
package main

import (
	"errors"
	"fmt"
	"log"
	"math"

	"github.com/rickykimani/zfactor/linalglite"
)

func main() {
	solveOnce()
	fmt.Println()

	reuseTheFactorisation()
	fmt.Println()

	singularSystem()
	fmt.Println()

	newtonIteration()
}

// solveOnce is the short form, for a system solved a single time.
func solveOnce() {
	// A material balance: three streams, three constraints.
	a := linalglite.NewFrom(3, 3, []float64{
		1, 1, 1,
		2, -1, 0,
		0, 3, -1,
	})
	b := []float64{100, 0, 0}

	x, err := linalglite.Solve(a, b)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Solving a 3x3 system once")
	fmt.Printf("  x = %.6f  %.6f  %.6f\n", x[0], x[1], x[2])

	// The useful check on a solution is whether it satisfies the
	// equations, not how it was obtained.
	residual, err := linalglite.Residual(a, x, b)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("  largest |A·x - b| = %.3e\n", residual)

	det, err := linalglite.Det(a)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("  determinant = %.6f\n", det)
	fmt.Println("\n  Orders up to three use closed forms rather than elimination,")
	fmt.Println("  which is several times faster and allocates once.")
}

// reuseTheFactorisation shows the form to prefer when one matrix is used
// against several right-hand sides.
//
// Factorising is the cubic part of the work and solving the quadratic
// part, so the saving grows with the number of right-hand sides.
func reuseTheFactorisation() {
	a := linalglite.NewFrom(4, 4, []float64{
		4, 1, 0, 0,
		1, 4, 1, 0,
		0, 1, 4, 1,
		0, 0, 1, 4,
	})

	lu, err := linalglite.Factorize(a)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("One factorisation, several right-hand sides")

	// SolveInto writes through a destination the caller owns, so nothing
	// is allocated per solve.
	x := make([]float64, 4)

	for _, b := range [][]float64{
		{1, 0, 0, 0},
		{0, 1, 0, 0},
		{1, 2, 3, 4},
	} {
		if err := lu.SolveInto(x, b); err != nil {
			log.Fatal(err)
		}

		fmt.Printf("  b = %v  ->  x = %8.5f %8.5f %8.5f %8.5f\n", b, x[0], x[1], x[2], x[3])
	}

	fmt.Printf("\n  determinant from the factorisation = %.6f\n", lu.Det())
}

// singularSystem shows what a system without a unique solution returns.
//
// A vanishing pivot could be divided by to produce infinities that look
// like numbers. Reporting it instead means a caller cannot mistake the
// result for an answer.
func singularSystem() {
	// The third row is the sum of the first two, so the rows are
	// dependent.
	a := linalglite.NewFrom(3, 3, []float64{
		1, 2, 3,
		4, 5, 6,
		5, 7, 9,
	})
	b := []float64{1, 1, 1}

	fmt.Println("A singular system")

	_, err := linalglite.Solve(a, b)

	if errors.Is(err, linalglite.ErrSingular) {
		fmt.Printf("  Solve reports: %v\n", err)
	}

	// A determinant, unlike a solution, exists for such a matrix and is
	// zero, so it is returned as a value.
	det, err := linalglite.Det(a)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Printf("  Det returns %.1f, with no error: a singular matrix has a determinant\n", det)
}

// newtonIteration is the pattern the package was written for.
//
// Each step evaluates a residual and a Jacobian, solves for a correction
// and applies it. The Jacobian changes every step, so Refactorize
// replaces the factorisation in the storage already held and the whole
// loop allocates nothing.
//
// The system solved here is a small nonlinear pair:
//
//	x² + y² = 25
//	x − y   = 1
//
// whose positive solution is x = 4, y = 3.
func newtonIteration() {
	const tolerance = 1e-12

	jacobian := linalglite.New(2, 2)
	residual := make([]float64, 2)
	step := make([]float64, 2)

	// A guess away from the answer, so the iteration has work to do.
	x, y := 1.0, 0.0

	evaluate := func() {
		residual[0] = x*x + y*y - 25
		residual[1] = x - y - 1

		row := jacobian.Row(0)
		row[0], row[1] = 2*x, 2*y

		row = jacobian.Row(1)
		row[0], row[1] = 1, -1
	}

	evaluate()

	lu, err := linalglite.Factorize(jacobian)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println("Newton's method on a nonlinear pair")
	fmt.Println("  step        x           y        |residual|")

	for iteration := range 12 {
		norm := math.Max(math.Abs(residual[0]), math.Abs(residual[1]))

		fmt.Printf("  %4d   %9.6f   %9.6f   %.3e\n", iteration, x, y, norm)

		if norm < tolerance {
			break
		}

		// The Jacobian has changed, so replace the factorisation without
		// allocating afresh.
		if err := lu.Refactorize(jacobian); err != nil {
			log.Fatal(err)
		}

		if err := lu.SolveInto(step, residual); err != nil {
			log.Fatal(err)
		}

		x -= step[0]
		y -= step[1]

		evaluate()
	}

	fmt.Printf("\n  converged to x = %.12f, y = %.12f (exactly 4 and 3)\n", x, y)
	fmt.Println("  Refactorize and SolveInto allocate nothing, so the loop")
	fmt.Println("  runs at a steady cost however many steps it takes.")
}
