package modified_raoult

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor/antoine"
)

// The methanol(1)/methyl acetate(2) system of Smith, Van Ness & Abbott,
// Example 13.1.
//
// The published Antoine correlations are expressed with T in kelvin,
//
//	ln P1_sat = 16.59158 - 3643.31/(T - 33.424)
//	ln P2_sat = 14.25326 - 2665.54/(T - 53.424)
//
// while this library evaluates the same form with t in °C. Since the
// temperature appears only as (T + C), the conversion is absorbed
// entirely into C: C_celsius = C_kelvin + 273.15. A and B are unchanged.
func methanolAntoine() *antoine.Antoine {
	return &antoine.Antoine{
		Name:    "methanol",
		A:       16.59158,
		B:       3643.31,
		C:       273.15 - 33.424,
		Range:   antoine.TempRange{Low: 0, High: 150},
		Formula: "CH4O",
	}
}

func methylAcetateAntoine() *antoine.Antoine {
	return &antoine.Antoine{
		Name:    "methyl acetate",
		A:       14.25326,
		B:       2665.54,
		C:       273.15 - 53.424,
		Range:   antoine.TempRange{Low: 0, High: 150},
		Formula: "C3H6O2",
	}
}

func exampleModels() []antoine.Model {
	return []antoine.Model{methanolAntoine(), methylAcetateAntoine()}
}

// The activity coefficients follow the symmetric two-parameter Margules
// model with a temperature-dependent parameter
//
//	A = 2.771 - 0.00523 T   (T in kelvin).
//
// The isothermal cases evaluate A once at the stated temperature.
func margulesA(tKelvin float64) float64 {
	return 2.771 - 0.00523*tKelvin
}

// TestAzeotropePExample13_1 checks the azeotropic pressure and
// composition of part (e): T = 318.15 K (45 °C).
//
// The published answers are x1 = y1 = 0.325 and P = 73.8 kPa.
func TestAzeotropePExample13_1(t *testing.T) {
	const (
		tKelvin  = 318.15
		tCelsius = tKelvin - 273.15

		wantX1 = 0.325
		wantP  = 73.8
	)

	a := margulesA(tKelvin)

	got, err := AzeotropeP(MixtureInput{
		T:        tCelsius,
		Antoine:  exampleModels(),
		Activity: Margules{A12: a, A21: a},
	})
	if err != nil {
		t.Fatalf("AzeotropeP returned an unexpected error: %v", err)
	}

	if len(got) != 1 {
		t.Fatalf("expected a single azeotrope; got %d", len(got))
	}

	res := got[0]

	if math.Abs(res.X[0]-wantX1) > 5e-3 {
		t.Errorf("azeotropic composition x1 = %.4f; want %.3f", res.X[0], wantX1)
	}

	if math.Abs(res.P-wantP) > 0.1 {
		t.Errorf("azeotropic pressure P = %.2f kPa; want %.1f kPa", res.P, wantP)
	}

	// At an azeotrope the phase compositions coincide.
	for i := range res.X {
		if math.Abs(res.X[i]-res.Y[i]) > 1e-9 {
			t.Errorf(
				"component %d: x = %.6f but y = %.6f; compositions must be equal",
				i, res.X[i], res.Y[i],
			)
		}
	}

	if res.T != tCelsius {
		t.Errorf("temperature = %.4f °C; want %.4f °C", res.T, tCelsius)
	}
}

// TestAzeotropeTRoundTrip checks the isobaric solver against the
// isothermal one: fixing the pressure at the azeotropic pressure found
// for 318.15 K must recover that temperature and composition.
func TestAzeotropeTRoundTrip(t *testing.T) {
	const (
		tKelvin  = 318.15
		tCelsius = tKelvin - 273.15
	)

	a := margulesA(tKelvin)

	isothermal, err := AzeotropeP(MixtureInput{
		T:        tCelsius,
		Antoine:  exampleModels(),
		Activity: Margules{A12: a, A21: a},
	})
	if err != nil {
		t.Fatalf("AzeotropeP returned an unexpected error: %v", err)
	}
	if len(isothermal) != 1 {
		t.Fatalf("expected a single azeotrope; got %d", len(isothermal))
	}

	isobaric, err := AzeotropeT(MixtureInput{
		P:        isothermal[0].P,
		Antoine:  exampleModels(),
		Activity: Margules{A12: a, A21: a},
	})
	if err != nil {
		t.Fatalf("AzeotropeT returned an unexpected error: %v", err)
	}
	if len(isobaric) != 1 {
		t.Fatalf("expected a single azeotrope; got %d", len(isobaric))
	}

	if math.Abs(isobaric[0].T-tCelsius) > 1e-2 {
		t.Errorf(
			"azeotropic temperature = %.4f °C; want %.4f °C",
			isobaric[0].T, tCelsius,
		)
	}

	if math.Abs(isobaric[0].X[0]-isothermal[0].X[0]) > 1e-3 {
		t.Errorf(
			"azeotropic composition x1 = %.4f; want %.4f",
			isobaric[0].X[0], isothermal[0].X[0],
		)
	}
}

// TestAzeotropePDouble checks that both azeotropes of a double
// azeotropic system are found.
//
// This is a regression test for the sweep introduced in
// azeotropicCompositions. With an even number of roots the residual has
// the same sign at both pure-component limits, so comparing only the
// endpoints reports no azeotrope at all.
//
// The system is constructed rather than measured: with the asymmetric
// Margules parameters below, the azeotropic condition reduces to
//
//	1 + 2x₁ - 6x₁² = ln(P₂_sat/P₁_sat)
//
// whose two roots lie inside (0, 1).
func TestAzeotropePDouble(t *testing.T) {
	const (
		wantFirst  = 0.046557
		wantSecond = 0.286776
	)

	got, err := AzeotropeP(SaturationPressureInput{
		T:        45,
		PSats:    []float64{30.0, 88.35},
		Activity: Margules{A12: 1.0, A21: 3.0},
	})
	if err != nil {
		t.Fatalf("AzeotropeP returned an unexpected error: %v", err)
	}

	if len(got) != 2 {
		t.Fatalf("expected two azeotropes; got %d", len(got))
	}

	if math.Abs(got[0].X[0]-wantFirst) > 1e-3 {
		t.Errorf("first azeotrope x1 = %.6f; want %.6f", got[0].X[0], wantFirst)
	}

	if math.Abs(got[1].X[0]-wantSecond) > 1e-3 {
		t.Errorf("second azeotrope x1 = %.6f; want %.6f", got[1].X[0], wantSecond)
	}

	if got[0].X[0] >= got[1].X[0] {
		t.Errorf("azeotropes must be ordered by increasing x1: got %v", got)
	}

	// Each azeotrope must satisfy P = γi Pi_sat for both components.
	for i, res := range got {
		gamma, err := activityCoefficients(
			Margules{A12: 1.0, A21: 3.0}, res.T, res.X,
		)
		if err != nil {
			t.Fatalf("activity coefficients failed: %v", err)
		}

		p1 := gamma[0] * 30.0
		p2 := gamma[1] * 88.35

		if math.Abs(p1-p2) > 1e-3 {
			t.Errorf(
				"azeotrope %d at x1 = %.6f: γ1P1sat = %.4f but γ2P2sat = %.4f",
				i, res.X[0], p1, p2,
			)
		}
	}
}

// TestAzeotropeNoAzeotrope checks that an ideal solution, whose
// relative volatility never reaches unity, is reported as having no
// azeotrope rather than as a solver failure.
func TestAzeotropeNoAzeotrope(t *testing.T) {
	_, err := AzeotropeP(MixtureInput{
		T:       45,
		Antoine: exampleModels(),
		// A = 0 gives unit activity coefficients (Raoult's law).
		Activity: Margules{A12: 0, A21: 0},
	})

	if !errors.Is(err, ErrNoAzeotrope) {
		t.Errorf("expected ErrNoAzeotrope for an ideal solution; got %v", err)
	}
}

// TestAzeotropeRejectsNonBinary checks that mixtures other than binaries
// are refused rather than silently mishandled.
func TestAzeotropeRejectsNonBinary(t *testing.T) {
	a := margulesA(318.15)

	_, err := AzeotropeP(SaturationPressureInput{
		T:        45,
		PSats:    []float64{44.51, 65.64, 50.0},
		Activity: Margules{A12: a, A21: a},
	})

	if err == nil {
		t.Fatal("expected an error for a ternary mixture; got nil")
	}

	if errors.Is(err, ErrNoAzeotrope) {
		t.Errorf("ternary input should be rejected as invalid, not as ErrNoAzeotrope")
	}
}
