// Package linalglite provides the dense linear algebra that engineering
// calculations keep needing: factorising a small matrix and solving
// against it, repeatedly, inside an iteration.
//
// It is deliberately narrow. A Newton step on a multicomponent
// equilibrium, a least-squares fit of model parameters to measured data,
// the extent of several simultaneous reactions — these want a few
// operations on matrices of order two to a few hundred, called often
// enough that allocation in the loop matters. Anything wider belongs in
// a full library such as gonum.
//
// # Solving a system
//
// The short form solves once and discards the work:
//
//	x, err := linalglite.Solve(a, b)
//
// Iterations should keep the factorisation instead, since the same matrix
// is often used against several right-hand sides:
//
//	lu, err := linalglite.Factorize(a)
//	if err != nil {
//	    return err
//	}
//
//	for range steps {
//	    if err := lu.SolveInto(x, b); err != nil {
//	        return err
//	    }
//	    // update b from x ...
//	}
//
// SolveInto writes through a destination the caller owns and allocates
// nothing. Where the matrix itself changes each step, as in a Newton
// iteration on a Jacobian, Refactorize reuses the storage already held:
//
//	for range steps {
//	    fillJacobian(a)
//
//	    if err := lu.Refactorize(a); err != nil {
//	        return err
//	    }
//
//	    if err := lu.SolveInto(step, residual); err != nil {
//	        return err
//	    }
//	}
//
// # Storage
//
// A Dense matrix holds its elements in a single slice, in row-major
// order, so that walking a row walks memory in order. Matrices are
// small enough here that one contiguous allocation and predictable
// access matter more than anything else.
//
// # Accuracy
//
// Factorisation uses partial pivoting: each column is searched for its
// largest remaining entry, which is moved to the diagonal before the
// column is eliminated. Without it, a small pivot divides the rest of
// the matrix and loses precision that no later step recovers. Matrices
// of order one, two and three are solved by their closed forms instead,
// which are both faster and exact to rounding.
//
// A matrix whose pivot vanishes is reported as singular rather than
// producing an infinity, since a system with no unique solution should
// not be answered with one.
package linalglite
