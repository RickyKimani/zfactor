package cubic

import (
	"errors"
	"math"
)

// LogFugacity calculates the natural logarithm of the fugacity coefficient.
// Z is the compressibility factor (PV/RT).
// A and B are the dimensionless EOS parameters: A = aP/(RT)^2, B = bP/RT.
func LogFugacity(cfg *EOSCfg, Z, A, B float64) float64 {
	sigma := cfg.Type.Params().Sigma
	epsilon := cfg.Type.Params().Epsilon

	// Generic cubic EOS fugacity coefficient
	// ln(phi) = Z - 1 - ln(Z - B) + (A / (B * (epsilon - sigma))) * ln((Z + sigma*B) / (Z + epsilon*B))

	term1 := Z - 1 - math.Log(Z-B)

	var term2 float64
	diff := epsilon - sigma
	if math.Abs(diff) < 1e-9 {
		// Degenerate case (e.g. VdW)
		// For VdW, the integral term is -A/Z.
		term2 = -A / Z
	} else {
		term2 = (A / (B * diff)) * math.Log((Z+sigma*B)/(Z+epsilon*B))
	}

	return term1 + term2
}

// SaturationPressure calculates the saturation pressure at a given temperature T.
// It uses the Wilson equation for the initial guess and iterates using the equal fugacity condition.
func SaturationPressure(cfg *EOSCfg, T float64) (float64, error) {
	if T >= cfg.Tc {
		return cfg.Pc, nil
	}

	// Initial guess using Wilson equation
	Tr := T / cfg.Tc
	P := cfg.Pc * math.Exp(5.373*(1+cfg.Acentric)*(1-1/Tr))

	for range 100 {
		// Update cfg with new P
		iterCfg := *cfg
		iterCfg.P = P
		iterCfg.T = T

		// Solve for volume
		volRes, err := SolveForVolume(&iterCfg)
		if err != nil {
			return 0, err
		}

		roots := volRes.Clean()

		// If we don't have 3 roots, we are likely outside the two-phase region (P is too high or too low).
		// We need to adjust P to find the 3-root region.
		if len(roots) < 3 {
			if len(roots) == 0 {
				return 0, errors.New("no real roots found")
			}

			// Heuristic: Check compressibility Z
			// If Z is small (liquid-like), P is likely too high -> decrease P
			// If Z is large (vapor-like), P is likely too low -> increase P
			// However, simply adjusting P blindly is dangerous.
			// Let's try to nudge P towards the Wilson guess if we drifted, or just fail gracefully.
			// For now, let's try to be robust: if we are close to the solution, we should have 3 roots.
			// If we don't, maybe our step was too big.

			// Alternative: If we have 1 root, we can't calculate fugacity difference.
			// But maybe we can use the single root to estimate? No.

			// Let's try to reduce P if V is small, increase if V is large?
			V := roots[0]
			b := volRes.B
			if V < 2*b {
				// Liquid-like, P too high?
				P = P * 0.9
			} else {
				// Vapor-like, P too low?
				P = P * 1.1
			}
			continue
		}

		Vl := roots[0]
		Vv := roots[len(roots)-1]

		// Calculate dimensionless parameters
		RT := cfg.R * T
		Adim := volRes.A * P / (RT * RT)
		Bdim := volRes.B * P / RT

		Zl := P * Vl / RT
		Zv := P * Vv / RT

		// Protect against invalid log arguments
		if Zl <= Bdim || Zv <= Bdim {
			// Should not happen for valid roots > b
			// But if it does, perturb P
			P = P * 0.95
			continue
		}

		phil := LogFugacity(&iterCfg, Zl, Adim, Bdim)
		phiv := LogFugacity(&iterCfg, Zv, Adim, Bdim)

		// Check convergence
		if math.Abs(phil-phiv) < 1e-8 {
			return P, nil
		}

		// Update P
		// P_new = P_old * exp(ln_phi_l - ln_phi_v)
		// Dampen the update to avoid oscillations
		ratio := math.Exp(phil - phiv)

		// Limit the step size
		if ratio > 1.2 {
			ratio = 1.2
		} else if ratio < 0.8 {
			ratio = 0.8
		}

		P = P * ratio
	}

	return 0, errors.New("saturation pressure did not converge")
}
