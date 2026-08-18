package cubic

import (
	"errors"
	"fmt"
	"math"
)

// ErrCompressibilityTooSmall is returned when the compressibility factor does
// not exceed the dimensionless parameter B. The residual properties both take
// ln(Z - B), so a Z at or below B has no residual property rather than a
// large one, and reporting it keeps a caller from reading an infinity as a
// number.
var ErrCompressibilityTooSmall = errors.New("compressibility factor must exceed B")

// LnAlphaDeriver is implemented by an equation of state that knows the
// temperature derivative of its own alpha function.
//
// The residual properties need
//
//	d ln α(Tr) / d ln Tr
//
// and nothing else in the package does. Putting it on EOSType would break
// every equation of state implemented outside this package for the sake of
// two functions, so it is asked for separately: an equation that provides it
// is differentiated exactly, and one that does not is differentiated
// numerically. The four equations here all provide it.
type LnAlphaDeriver interface {
	// DLnAlphaDLnTr returns d ln α / d ln Tr at the given reduced
	// temperature and acentric factor.
	DLnAlphaDLnTr(tr, w float64) float64
}

// DLnAlphaDLnTr for van der Waals, whose alpha is 1 at every temperature and
// so has no temperature dependence to contribute.
func (*VdW) DLnAlphaDLnTr(tr, w float64) float64 {
	return 0
}

// DLnAlphaDLnTr for Redlich-Kwong, where α = Tr^(-1/2), so ln α is
// -(1/2) ln Tr and the derivative is a constant.
func (*RK) DLnAlphaDLnTr(tr, w float64) float64 {
	return -0.5
}

// DLnAlphaDLnTr for Soave-Redlich-Kwong.
func (*SRK) DLnAlphaDLnTr(tr, w float64) float64 {
	return soaveLnAlphaDeriv(srkM(w), tr)
}

// DLnAlphaDLnTr for Peng-Robinson, which shares Soave's form of alpha and
// differs only in the polynomial for m.
func (*PR) DLnAlphaDLnTr(tr, w float64) float64 {
	return soaveLnAlphaDeriv(prM(w), tr)
}

// soaveLnAlphaDeriv returns d ln α / d ln Tr for the Soave form of alpha,
//
//	α = u²  where  u = 1 + m(1 - √Tr)
//
// which SRK and PR share. Differentiating,
//
//	du/dTr        = -m / (2√Tr)
//	dα/dTr        = 2u · du/dTr = -m·u / √Tr
//	d ln α/d ln Tr = (Tr/α)(dα/dTr) = -m·√Tr / u
//
// At Tr = 1 this is -m, which is the value quoted for the slope of alpha at
// the critical temperature and the check the tests pin it against.
func soaveLnAlphaDeriv(m, tr float64) float64 {
	if tr <= 0 {
		return math.NaN()
	}

	root := math.Sqrt(tr)
	u := 1 + m*(1-root)

	// u is zero where alpha is, and ln α is undefined there, so its
	// derivative is reported as undefined rather than as an infinity.
	if u == 0 {
		return math.NaN()
	}

	return -m * root / u
}

// lnAlphaDeriv returns d ln α / d ln Tr for any equation of state, exactly
// when it can and by central difference otherwise.
//
// The difference is taken in ln Tr directly, which is the variable being
// differentiated against, so the step is relative and suits alpha functions
// that vary over orders of magnitude in temperature.
func lnAlphaDeriv(eos EOSType, tr, w float64) float64 {
	if exact, ok := eos.(LnAlphaDeriver); ok {
		return exact.DLnAlphaDLnTr(tr, w)
	}

	// A central difference converges as h², and its total error is least
	// near the cube root of the machine epsilon.
	const h = 1e-5

	up := eos.Alpha(tr*math.Exp(h), w)
	down := eos.Alpha(tr*math.Exp(-h), w)

	// ln alpha is undefined for a non-positive alpha, so there is nothing
	// to difference.
	if up <= 0 || down <= 0 {
		return math.NaN()
	}

	return (math.Log(up) - math.Log(down)) / (2 * h)
}

// qTimesI returns the product q·I that appears in both residual properties,
// where
//
//	q = Ψ α(Tr) / (Ω Tr)   and   I = (1/(σ-ε)) · ln[(Z + σB) / (Z + εB)]
//
// B is the dimensionless b P / (RT), so q is A/B for the dimensionless
// A = a P / (RT)².
//
// When σ and ε are equal the expression for I is a difference of equal
// logarithms over a vanishing denominator. Its limit is finite,
//
//	I → B / (Z + εB)
//
// which is the form van der Waals needs.
func qTimesI(cfg *EOSCfg, Z, A, B float64) (float64, error) {
	params := cfg.Type.Params()
	sigma, epsilon := params.Sigma, params.Epsilon

	if B <= 0 {
		return 0, fmt.Errorf("B must be greater than 0, got %g", B)
	}

	q := A / B

	// The gap between sigma and epsilon is a property of the equation, not
	// of the state, so this selects a branch rather than guarding a
	// near-singularity.
	gap := sigma - epsilon

	var I float64

	if gap == 0 {
		denom := Z + epsilon*B
		if denom == 0 {
			return 0, fmt.Errorf("%w: Z + εB is 0", ErrCompressibilityTooSmall)
		}

		I = B / denom
	} else {
		upper := Z + sigma*B
		lower := Z + epsilon*B

		// Either being non-positive puts the ratio outside the logarithm's
		// domain, which happens for a Z too small for the equation rather
		// than through any failure here.
		if upper <= 0 || lower <= 0 {
			return 0, fmt.Errorf("%w: Z + σB = %g and Z + εB = %g must both be positive",
				ErrCompressibilityTooSmall, upper, lower)
		}

		I = math.Log(upper/lower) / gap
	}

	return q * I, nil
}

// ResidualEnthalpy returns the residual enthalpy of a state, as the
// dimensionless group H^R/(R Tc).
//
// The underlying equation is written in terms of H^R/(RT),
//
//	H^R / (R T) = Z - 1 + [d ln α(Tr)/d ln Tr - 1] q I
//
// and the two differ by a factor of the reduced temperature,
//
//	H^R / (R Tc) = Tr · H^R / (R T)
//
// which is applied here so that this agrees with the correlations of the same
// name in the abbott and lee-kesler packages. All three return H^R/(R Tc),
// and all three return S^R/R for the entropy. Multiply by R Tc to recover an
// enthalpy in the units of R.
//
// Z is the compressibility factor of the phase in question, so a state with
// both a liquid and a vapour root has a residual enthalpy for each and the
// caller chooses by which Z it passes. A and B are the dimensionless
// parameters a P / (RT)² and b P / (RT), the same pair LogFugacity takes.
func ResidualEnthalpy(cfg *EOSCfg, Z, A, B float64) (float64, error) {
	if cfg == nil {
		return 0, errors.New("configuration error: config cannot be nil")
	}
	if cfg.Type == nil {
		return 0, errors.New("configuration error: 'Type' field (EOS model) is required")
	}
	if cfg.Tc <= 0 {
		return 0, errors.New("configuration error: critical temperature must be greater than 0")
	}

	qI, err := qTimesI(cfg, Z, A, B)
	if err != nil {
		return 0, err
	}

	tr := cfg.T / cfg.Tc

	overRT := Z - 1 + (lnAlphaDeriv(cfg.Type, tr, cfg.Acentric)-1)*qI

	return tr * overRT, nil
}

// ResidualEntropy returns the residual entropy of a state, as the
// dimensionless group
//
//	S^R / R = ln(Z - B) + [d ln α(Tr)/d ln Tr] q I
//
// The arguments carry the same meaning as in ResidualEnthalpy. Multiply by R
// to recover an entropy in the units of R.
func ResidualEntropy(cfg *EOSCfg, Z, A, B float64) (float64, error) {
	if cfg == nil {
		return 0, errors.New("configuration error: config cannot be nil")
	}
	if cfg.Type == nil {
		return 0, errors.New("configuration error: 'Type' field (EOS model) is required")
	}
	if cfg.Tc <= 0 {
		return 0, errors.New("configuration error: critical temperature must be greater than 0")
	}

	qI, err := qTimesI(cfg, Z, A, B)
	if err != nil {
		return 0, err
	}

	// The logarithm is what makes Z > B a requirement rather than a
	// preference: at Z = B the state is at the hard-sphere limit of the
	// equation, where the residual entropy diverges.
	if Z <= B {
		return 0, fmt.Errorf("%w: Z = %g, B = %g", ErrCompressibilityTooSmall, Z, B)
	}

	tr := cfg.T / cfg.Tc

	return math.Log(Z-B) + lnAlphaDeriv(cfg.Type, tr, cfg.Acentric)*qI, nil
}
