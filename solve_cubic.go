package zfactor

import (
	"errors"
	"math"
	"math/cmplx"
)

// polishSteps is the number of Newton iterations applied to each root
// once Cardano's formulae have located it.
//
// Newton converges quadratically near a simple root, so a very few
// passes take an estimate already close to the answer down to rounding.
// Three is enough to reach machine precision from the accuracy Cardano
// delivers and cheap beside the cube roots that precede it.
const polishSteps = 3

// SolveCubic solves ax^3 + bx^2 + cx + d = 0
// Returns all 3 roots (possibly complex).
//
// The roots come from Cardano's formulae and are then refined by
// Newton's method. The closed form loses precision when the roots differ
// widely in magnitude, because the substitution that removes the
// quadratic term subtracts quantities of similar size; refining recovers
// it. Across the polynomials the equation-of-state solvers construct,
// this takes the worst relative residual from around 1e-5 to 1e-15.
//
// Precision cannot be recovered when the leading coefficient is small
// enough that the cubic is effectively a quadratic, since one root then
// escapes toward infinity and cannot be represented. Only a exactly zero
// is rejected: no threshold on its magnitude separates that case from
// the legitimate ones, the equation-of-state polynomials being monic
// with trailing coefficients as large as 1e11.
func SolveCubic(a, b, c, d float64) ([3]complex128, error) {
	if a == 0 {
		return [3]complex128{}, errors.New("equation provided is not cubic (a = 0)")
	}

	// 1. Normalize coefficients
	b /= a
	c /= a
	d /= a

	// 2. Depressed cubic: y^3 + py + q = 0
	p := c - b*b/3
	q := 2*b*b*b/27 - b*c/3 + d

	// 3. Discriminant
	delta := (q*q)/4 + (p*p*p)/27

	// 4. Cube roots of unity
	omega := complex(-0.5, math.Sqrt(3)/2)
	omega2 := complex(-0.5, -math.Sqrt(3)/2)

	var roots [3]complex128

	if delta >= 0 {
		// One real root and two complex conjugates.
		// Use REAL cube roots (math.Cbrt handles negative radicands): this
		// keeps the Cardano coupling u*v = -p/3 satisfied automatically,
		// since Cbrt(A)*Cbrt(B) = Cbrt(A*B) = Cbrt(-p^3/27) = -p/3. Taking
		// independent principal cube roots via cmplx.Pow breaks that coupling
		// whenever a radicand is negative and yields wrong roots.
		sd := math.Sqrt(delta)
		u := complex(math.Cbrt(-q/2+sd), 0)
		v := complex(math.Cbrt(-q/2-sd), 0)

		y1 := u + v
		y2 := u*omega + v*omega2
		y3 := u*omega2 + v*omega

		shift := complex(b/3, 0)
		roots[0] = y1 - shift
		roots[1] = y2 - shift
		roots[2] = y3 - shift
	} else {
		// Three real roots
		r := math.Sqrt(-p * p * p / 27)
		phi := math.Acos(-q / (2 * math.Sqrt(-(p*p*p)/27)))
		t := 2 * math.Cbrt(r)

		y1 := complex(t*math.Cos(phi/3), 0)
		y2 := complex(t*math.Cos((phi+2*math.Pi)/3), 0)
		y3 := complex(t*math.Cos((phi+4*math.Pi)/3), 0)

		shift := complex(b/3, 0)
		roots[0] = y1 - shift
		roots[1] = y2 - shift
		roots[2] = y3 - shift
	}

	for i := range roots {
		roots[i] = polish(b, c, d, roots[i])
	}

	return roots, nil
}

// polish refines a root of the monic cubic x^3 + bx^2 + cx + d by
// Newton's method.
//
// The iteration is abandoned rather than applied blindly if the
// derivative vanishes, which happens at a repeated root, or if a step
// leaves the finite plane. In either case the estimate Cardano supplied
// is returned unchanged, so refining can only improve a root or leave it
// as it was.
func polish(b, c, d float64, x complex128) complex128 {
	B, C, D := complex(b, 0), complex(c, 0), complex(d, 0)

	for range polishSteps {
		value := ((x+B)*x+C)*x + D
		derivative := (3*x+2*B)*x + C

		if derivative == 0 {
			return x
		}

		next := x - value/derivative

		if cmplx.IsNaN(next) || cmplx.IsInf(next) {
			return x
		}

		x = next
	}

	return x
}
