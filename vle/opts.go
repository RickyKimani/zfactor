// Package vle provides shared types and utilities used by vapor-liquid
// equilibrium (VLE) calculations.
//
// The package contains common configuration shared by ideal and
// non-ideal VLE models, including numerical solver options. Specific
// equilibrium models are implemented by subpackages such as raoult and
// modified_raoult.
package vle

const (
	defaultRootTolerance = 1e-6
	defaultMaxIterations = 100
)

// SolverOptions configures the numerical algorithms used by iterative
// VLE calculations.
//
// A zero-valued SolverOptions is valid and causes the default
// convergence tolerance and maximum iteration count to be used.
type SolverOptions struct {
	// Tolerance specifies the convergence criterion used by iterative
	// solvers. If zero or negative, a default tolerance of 1×10⁻⁶ is
	// used.
	Tolerance float64

	// MaxIterations specifies the maximum number of iterations allowed
	// before a solver reports failure. If zero or negative, a default
	// limit of 10 iterations is used.
	MaxIterations int
}

// Tol returns the convergence tolerance, substituting the package
// default when no positive tolerance has been specified.
func (s SolverOptions) Tol() float64 {
	if s.Tolerance <= 0 {
		return defaultRootTolerance
	}
	return s.Tolerance
}

// MaxIter returns the maximum iteration count, substituting the
// package default when no positive limit has been specified.
func (s SolverOptions) MaxIter() int {
	if s.MaxIterations <= 0 {
		return defaultMaxIterations
	}
	return s.MaxIterations
}
