package zfactor

import (
	"errors"
	"math"
)

// SolveCubic solves ax^3 + bx^2 + cx + d = 0
// Returns all 3 roots (possibly complex).
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

	return roots, nil
}
