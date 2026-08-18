package cubic_test

import (
	"errors"
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/substance"
)

// R in bar*cm^3/(mol*K), the convention used across the package.
const residualR = 10 * zfactor.RSI

// eosState is a state reduced to what the residual properties need: the
// compressibility factor of one phase and the dimensionless parameters.
type eosState struct {
	cfg *cubic.EOSCfg
	Z   float64
	A   float64
	B   float64
}

// vaporState solves the equation of state at (T, P) and returns the largest
// real root, which is the vapour where more than one exists.
//
// It reproduces the dimensionless grouping the package works in: A = aP/(RT)^2
// and B = bP/RT, with Z = PV/RT.
func vaporState(t *testing.T, base *cubic.EOSCfg, T, P float64) eosState {
	t.Helper()

	cfg := *base
	cfg.T, cfg.P = T, P

	res, err := cubic.SolveForVolume(&cfg)
	if err != nil {
		t.Fatalf("SolveForVolume at T=%g P=%g: %v", T, P, err)
	}

	roots := res.Clean()
	if len(roots) == 0 {
		t.Fatalf("no real roots at T=%g P=%g", T, P)
	}

	RT := cfg.R * T

	return eosState{
		cfg: &cfg,
		Z:   P * roots[len(roots)-1] / RT,
		A:   res.A * P / (RT * RT),
		B:   res.B * P / RT,
	}
}

// equations returns the four equations of state against a substance whose
// acentric factor is not zero, so the terms carrying omega are exercised.
func equations(t *testing.T) map[string]*cubic.EOSCfg {
	t.Helper()

	sub := substance.Ethane
	Tc, Pc, w := sub.Critical.Tc, sub.Critical.Pc, sub.Acentric

	return map[string]*cubic.EOSCfg{
		"vdW": cubic.NewvdWCfg(0, 0, Tc, Pc, residualR),
		"RK":  cubic.NewRKCfg(0, 0, Tc, Pc, residualR),
		"SRK": cubic.NewSRKCfg(0, 0, Tc, Pc, w, residualR),
		"PR":  cubic.NewPRCfg(0, 0, Tc, Pc, w, residualR),
	}
}

// The residual Gibbs energy can be reached two ways, and they must agree:
//
//	G^R/RT = ln phi        and        G^R/RT = H^R/RT - S^R/R
//
// The identity is exact, and it holds whatever the alpha derivative is,
// because that term enters the enthalpy and the entropy with coefficients
// that differ by exactly one and so cancels. That makes this a check on the
// rest of both equations against LogFugacity, which computes its integral by
// a separate route.
func TestResidualPropertiesAgreeWithFugacity(t *testing.T) {
	states := []struct {
		name string
		T, P float64
	}{
		{"gas well above Tc", 400, 10},
		{"gas near Tc", 320, 30},
		{"compressed", 350, 80},
		{"low pressure", 300, 1},
	}

	for name, cfg := range equations(t) {
		for _, s := range states {
			t.Run(name+" "+s.name, func(t *testing.T) {
				state := vaporState(t, cfg, s.T, s.P)

				hr, err := cubic.ResidualEnthalpy(state.cfg, state.Z, state.A, state.B)
				if err != nil {
					t.Fatalf("ResidualEnthalpy: %v", err)
				}

				sr, err := cubic.ResidualEntropy(state.cfg, state.Z, state.A, state.B)
				if err != nil {
					t.Fatalf("ResidualEntropy: %v", err)
				}

				lnPhi := cubic.LogFugacity(state.cfg, state.Z, state.A, state.B)

				// ResidualEnthalpy is normalised by R Tc, and the identity
				// is in terms of H^R/RT, so the reduced temperature comes
				// back out here.
				tr := state.cfg.T / state.cfg.Tc
				overRT := hr / tr

				const tolerance = 1e-12

				if diff := math.Abs((overRT - sr) - lnPhi); diff > tolerance {
					t.Errorf("H^R/RT - S^R/R = %.15g, ln phi = %.15g, differ by %.3g",
						overRT-sr, lnPhi, diff)
				}
			})
		}
	}
}

// lnPhiAt returns ln phi of the vapour at a state, recomputing the
// temperature-dependent parameters of the equation from scratch.
func lnPhiAt(t *testing.T, base *cubic.EOSCfg, T, P float64) float64 {
	t.Helper()

	state := vaporState(t, base, T, P)

	return cubic.LogFugacity(state.cfg, state.Z, state.A, state.B)
}

// The alpha derivative cancels out of the identity above, so it needs a check
// of its own. Gibbs-Helmholtz supplies one:
//
//	(d(G^R/RT)/dT)_P = -H^R/(RT^2)   so   H^R/RT = -T (d ln phi/dT)_P
//
// Differentiating ln phi numerically against temperature therefore reproduces
// the residual enthalpy without using equation 13.75 at all, and the
// temperature dependence it picks up is exactly the alpha derivative, along
// with the way Z shifts as the equation is re-solved at each temperature.
func TestResidualEnthalpyMatchesGibbsHelmholtz(t *testing.T) {
	// A single-phase gas, so the root followed stays the same one as the
	// temperature is stepped and ln phi is smooth across the interval.
	const (
		T = 400.0
		P = 10.0
	)

	for name, cfg := range equations(t) {
		t.Run(name, func(t *testing.T) {
			state := vaporState(t, cfg, T, P)

			direct, err := cubic.ResidualEnthalpy(state.cfg, state.Z, state.A, state.B)
			if err != nil {
				t.Fatalf("ResidualEnthalpy: %v", err)
			}

			// A relative step, since the derivative is against temperature.
			h := T * 1e-5

			slope := (lnPhiAt(t, cfg, T+h, P) - lnPhiAt(t, cfg, T-h, P)) / (2 * h)
			viaGibbsHelmholtz := -T * slope

			// Gibbs-Helmholtz gives H^R/RT while ResidualEnthalpy returns
			// H^R/(R Tc), so compare on the former.
			tr := state.cfg.T / state.cfg.Tc
			overRT := direct / tr

			// The tolerance is set by the central difference, not by the
			// equations: its truncation error falls as h^2.
			const tolerance = 1e-6

			if diff := math.Abs(overRT - viaGibbsHelmholtz); diff > tolerance {
				t.Errorf("H^R/RT = %.12g from eq 13.75, %.12g from Gibbs-Helmholtz, differ by %.3g",
					overRT, viaGibbsHelmholtz, diff)
			}
		})
	}
}

// The alpha derivative is written in closed form for each equation. Comparing
// it against a numerical derivative of Alpha checks the differentiation
// itself, independently of any thermodynamics, and would catch a closed form
// that no longer matches the alpha it belongs to.
func TestLnAlphaDerivativeMatchesNumericalDifferentiation(t *testing.T) {
	sub := substance.Ethane
	Tc, Pc := sub.Critical.Tc, sub.Critical.Pc

	acentrics := []float64{0, 0.1, 0.3, 0.5}
	reduced := []float64{0.5, 0.7, 0.9, 1.0, 1.3, 2.0, 5.0}

	for _, w := range acentrics {
		eqs := map[string]cubic.EOSType{
			"vdW": cubic.NewvdWCfg(0, 0, Tc, Pc, residualR).Type,
			"RK":  cubic.NewRKCfg(0, 0, Tc, Pc, residualR).Type,
			"SRK": cubic.NewSRKCfg(0, 0, Tc, Pc, w, residualR).Type,
			"PR":  cubic.NewPRCfg(0, 0, Tc, Pc, w, residualR).Type,
		}

		for name, eos := range eqs {
			deriver, ok := eos.(cubic.LnAlphaDeriver)
			if !ok {
				t.Fatalf("%s does not provide its alpha derivative", name)
			}

			for _, tr := range reduced {
				exact := deriver.DLnAlphaDLnTr(tr, w)

				// Central difference in ln Tr, the variable being
				// differentiated against.
				const h = 1e-6

				up := eos.Alpha(tr*math.Exp(h), w)
				down := eos.Alpha(tr*math.Exp(-h), w)
				numerical := (math.Log(up) - math.Log(down)) / (2 * h)

				// Relative, since the derivative spans a wide range across
				// this grid.
				const tolerance = 1e-6

				if diff := math.Abs(exact - numerical); diff > tolerance*math.Max(1, math.Abs(exact)) {
					t.Errorf("%s w=%g Tr=%g: closed form %.12g, numerical %.12g, differ by %.3g",
						name, w, tr, exact, numerical, diff)
				}
			}
		}
	}
}

// At the critical temperature the Soave form reduces to -m, and the constant
// alphas reduce to values that can be written down. The polynomials for m are
// transcribed here rather than read from the package, so this is an
// independent statement of what they should be.
func TestLnAlphaDerivativeAtCriticalTemperature(t *testing.T) {
	sub := substance.Ethane
	Tc, Pc := sub.Critical.Tc, sub.Critical.Pc

	for _, w := range []float64{0, 0.1, 0.3} {
		cases := []struct {
			name string
			eos  cubic.EOSType
			want float64
		}{
			{"vdW", cubic.NewvdWCfg(0, 0, Tc, Pc, residualR).Type, 0},
			{"RK", cubic.NewRKCfg(0, 0, Tc, Pc, residualR).Type, -0.5},
			{
				"SRK",
				cubic.NewSRKCfg(0, 0, Tc, Pc, w, residualR).Type,
				-(0.480 + 1.574*w - 0.176*w*w),
			},
			{
				"PR",
				cubic.NewPRCfg(0, 0, Tc, Pc, w, residualR).Type,
				-(0.37464 + 1.54226*w - 0.26992*w*w),
			},
		}

		for _, c := range cases {
			deriver, ok := c.eos.(cubic.LnAlphaDeriver)
			if !ok {
				t.Fatalf("%s does not provide its alpha derivative", c.name)
			}

			got := deriver.DLnAlphaDLnTr(1.0, w)

			if math.Abs(got-c.want) > 1e-12 {
				t.Errorf("%s w=%g at Tr=1: got %.12g, want %.12g", c.name, w, got, c.want)
			}
		}
	}
}

// Both residual properties measure the departure of a real fluid from an
// ideal gas, so both must vanish as the pressure does. This limit is fixed by
// thermodynamics, not by any choice made here.
func TestResidualPropertiesVanishAsPressureFalls(t *testing.T) {
	const T = 400.0

	for name, cfg := range equations(t) {
		t.Run(name, func(t *testing.T) {
			var previousH, previousS float64

			for i, P := range []float64{1, 0.1, 0.01, 0.001} {
				state := vaporState(t, cfg, T, P)

				hr, err := cubic.ResidualEnthalpy(state.cfg, state.Z, state.A, state.B)
				if err != nil {
					t.Fatalf("ResidualEnthalpy at P=%g: %v", P, err)
				}

				sr, err := cubic.ResidualEntropy(state.cfg, state.Z, state.A, state.B)
				if err != nil {
					t.Fatalf("ResidualEntropy at P=%g: %v", P, err)
				}

				// Each tenfold drop in pressure should shrink both, since the
				// leading departure is linear in pressure.
				if i > 0 {
					if math.Abs(hr) >= math.Abs(previousH) {
						t.Errorf("H^R/RT did not shrink with pressure: %.6g at P=%g after %.6g",
							hr, P, previousH)
					}
					if math.Abs(sr) >= math.Abs(previousS) {
						t.Errorf("S^R/R did not shrink with pressure: %.6g at P=%g after %.6g",
							sr, P, previousS)
					}
				}

				previousH, previousS = hr, sr
			}

			// Shrinking is the weaker half of the statement. The leading
			// departure is linear in pressure, so the ratio of either
			// residual property to the pressure approaches a constant as
			// the pressure falls — checking that the ratio stops moving
			// pins the approach to the ideal gas, where an absolute
			// threshold would only say the numbers are small.
			var previousRatio float64

			for i, P := range []float64{0.01, 1e-3, 1e-4} {
				state := vaporState(t, cfg, T, P)

				hr, err := cubic.ResidualEnthalpy(state.cfg, state.Z, state.A, state.B)
				if err != nil {
					t.Fatalf("ResidualEnthalpy at P=%g: %v", P, err)
				}

				ratio := hr / P

				if i > 0 {
					// A tenfold drop in pressure should leave the ratio
					// where it was.
					const tolerance = 1e-3

					if change := math.Abs(ratio/previousRatio - 1); change > tolerance {
						t.Errorf("H^R/(R Tc P) moved from %.8g to %.8g as P fell to %g, a change of %.2g",
							previousRatio, ratio, P, change)
					}
				}

				previousRatio = ratio
			}
		})
	}
}

// A compressibility factor at or below B leaves ln(Z - B) undefined. The state
// is outside what the equation describes rather than the calculation having
// failed, so it is reported instead of returning an infinity that would read
// as a number.
func TestResidualPropertiesRejectCompressibilityBelowB(t *testing.T) {
	sub := substance.Ethane
	cfg := cubic.NewPRCfg(400, 10, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	const B = 0.05

	for _, Z := range []float64{B, B / 2, 0} {
		if _, err := cubic.ResidualEntropy(cfg, Z, 0.1, B); !errors.Is(err, cubic.ErrCompressibilityTooSmall) {
			t.Errorf("ResidualEntropy with Z=%g, B=%g: got %v, want ErrCompressibilityTooSmall", Z, B, err)
		}
	}
}

func TestResidualPropertiesRejectBadConfig(t *testing.T) {
	sub := substance.Ethane
	valid := cubic.NewPRCfg(400, 10, sub.Critical.Tc, sub.Critical.Pc, sub.Acentric, residualR)

	noType := *valid
	noType.Type = nil

	noTc := *valid
	noTc.Tc = 0

	cases := []struct {
		name string
		cfg  *cubic.EOSCfg
	}{
		{"nil config", nil},
		{"no equation of state", &noType},
		{"no critical temperature", &noTc},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, err := cubic.ResidualEnthalpy(c.cfg, 0.9, 0.1, 0.05); err == nil {
				t.Error("ResidualEnthalpy returned no error")
			}
			if _, err := cubic.ResidualEntropy(c.cfg, 0.9, 0.1, 0.05); err == nil {
				t.Error("ResidualEntropy returned no error")
			}
		})
	}
}

// Example 13.5 and Table 13.7 of Smith, Van Ness and Abbott: the residual
// enthalpy and entropy of n-butane gas at 500 K and 50 bar, worked for the
// Redlich-Kwong equation in the text and tabulated for all four.
//
// This is the only check here against published numbers rather than against
// an identity, so it is what confirms the convention as well as the algebra:
// an error repeated consistently in both equations would satisfy every
// identity in this file and still fail here.
//
// The book rounds its intermediates — Tr to 1.176, Pr to 1.317, Z to four
// figures — and quotes R as 8.314, so the comparison is relative rather than
// exact. Agreement is well inside a tenth of a percent on all twelve values.
func TestResidualPropertiesAgainstTable13_7(t *testing.T) {
	// The book's gas constant, so the comparison does not turn on a
	// difference in R.
	const (
		bookR = 8.314 // J/(mol*K)
		T     = 500.0 // K
		P     = 50.0  // bar
	)

	sub := substance.NButane
	Tc, Pc, w := sub.Critical.Tc, sub.Critical.Pc, sub.Acentric

	cases := []struct {
		name string
		cfg  *cubic.EOSCfg
		// Published in Table 13.7.
		wantZ float64
		wantH float64 // J/mol
		wantS float64 // J/(mol*K)
	}{
		{"vdW", cubic.NewvdWCfg(T, P, Tc, Pc, residualR), 0.6608, -3937, -5.424},
		{"RK", cubic.NewRKCfg(T, P, Tc, Pc, residualR), 0.6850, -4505, -6.546},
		{"SRK", cubic.NewSRKCfg(T, P, Tc, Pc, w, residualR), 0.7222, -4824, -7.413},
		{"PR", cubic.NewPRCfg(T, P, Tc, Pc, w, residualR), 0.6907, -4988, -7.426},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			state := vaporState(t, c.cfg, T, P)

			hr, err := cubic.ResidualEnthalpy(state.cfg, state.Z, state.A, state.B)
			if err != nil {
				t.Fatalf("ResidualEnthalpy: %v", err)
			}

			sr, err := cubic.ResidualEntropy(state.cfg, state.Z, state.A, state.B)
			if err != nil {
				t.Fatalf("ResidualEntropy: %v", err)
			}

			// H^R/(R Tc) back to J/mol, and S^R/R to J/(mol*K).
			gotH := hr * bookR * Tc
			gotS := sr * bookR

			// A tenth of a percent, which the rounding in the published
			// intermediates accounts for several times over.
			const tolerance = 1e-3

			checks := []struct {
				label     string
				got, want float64
			}{
				{"Z", state.Z, c.wantZ},
				{"H^R (J/mol)", gotH, c.wantH},
				{"S^R (J/mol/K)", gotS, c.wantS},
			}

			for _, check := range checks {
				if relative := math.Abs(check.got/check.want - 1); relative > tolerance {
					t.Errorf("%s: got %.6g, book gives %.6g, off by %.3g%%",
						check.label, check.got, check.want, 100*relative)
				}
			}
		})
	}
}

// The worked solution in the text publishes its intermediates for the
// Redlich-Kwong case, which pins the pieces of the calculation rather than
// only the result. Each is quoted to the digits shown there.
func TestExample13_5Intermediates(t *testing.T) {
	const (
		T = 500.0
		P = 50.0
	)

	sub := substance.NButane
	cfg := cubic.NewRKCfg(T, P, sub.Critical.Tc, sub.Critical.Pc, residualR)

	state := vaporState(t, cfg, T, P)

	// beta = Omega Pr / Tr, which is the dimensionless B the package works in.
	const wantBeta = 0.09703

	if relative := math.Abs(state.B/wantBeta - 1); relative > 1e-3 {
		t.Errorf("beta: got %.6g, book gives %.6g, off by %.3g%%", state.B, wantBeta, 100*relative)
	}

	// q = Psi alpha / (Omega Tr), which is A/B.
	const wantQ = 3.8689

	if q := state.A / state.B; math.Abs(q/wantQ-1) > 1e-3 {
		t.Errorf("q: got %.6g, book gives %.6g, off by %.3g%%", q, wantQ, 100*math.Abs(q/wantQ-1))
	}

	// I = ln[(Z + sigma beta)/(Z + epsilon beta)] / (sigma - epsilon), which
	// for Redlich-Kwong is ln[(Z + beta)/Z].
	const wantI = 0.13247

	gotI := math.Log((state.Z + state.B) / state.Z)

	if relative := math.Abs(gotI/wantI - 1); relative > 1e-3 {
		t.Errorf("I: got %.6g, book gives %.6g, off by %.3g%%", gotI, wantI, 100*relative)
	}

	// And the dimensionless groups the text prints before converting them.
	const (
		wantHOverRT = -1.0838
		wantSOverR  = -0.78735
	)

	hr, err := cubic.ResidualEnthalpy(state.cfg, state.Z, state.A, state.B)
	if err != nil {
		t.Fatalf("ResidualEnthalpy: %v", err)
	}

	// Undo the R Tc normalisation to compare against the printed H^R/RT.
	tr := T / sub.Critical.Tc
	gotHOverRT := hr / tr

	if relative := math.Abs(gotHOverRT/wantHOverRT - 1); relative > 1e-3 {
		t.Errorf("H^R/RT: got %.6g, book gives %.6g, off by %.3g%%",
			gotHOverRT, wantHOverRT, 100*relative)
	}

	sr, err := cubic.ResidualEntropy(state.cfg, state.Z, state.A, state.B)
	if err != nil {
		t.Fatalf("ResidualEntropy: %v", err)
	}

	if relative := math.Abs(sr/wantSOverR - 1); relative > 1e-3 {
		t.Errorf("S^R/R: got %.6g, book gives %.6g, off by %.3g%%", sr, wantSOverR, 100*relative)
	}
}
