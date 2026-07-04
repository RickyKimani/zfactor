package modified_raoult

import (
	"github.com/rickykimani/zfactor/activity"
	"github.com/rickykimani/zfactor/activity/margules"
	"github.com/rickykimani/zfactor/activity/nrtl"
	vanlaar "github.com/rickykimani/zfactor/activity/van-laar"
	"github.com/rickykimani/zfactor/activity/wilson"
)

// ActivityModel constructs an activity-coefficient model for use by the
// modified Raoult's law VLE calculations.
//
// The wrappers provided by this package allow VLE calculations to remain
// independent of the implementations in the activity package while
// exposing a simple API suited for vapor-liquid equilibrium calculations.
type ActivityModel interface {
	// Model constructs the underlying activity model.
	Model() activity.Model
}

// Wilson supplies the parameters required to construct a Wilson
// activity-coefficient model.
type Wilson struct {
	// V contains the pure-component molar liquid volumes.
	V []float64

	// Interaction contains the Wilson interaction energy matrix.
	// The values are typically expressed in J/mol.
	Interaction [][]float64
}

// Model constructs the corresponding Wilson activity model.
func (w Wilson) Model() activity.Model {
	return wilson.Wilson{
		V:           w.V,
		Interaction: w.Interaction,
	}
}

// NRTL supplies the parameters required to construct an NRTL
// activity-coefficient model.
type NRTL struct {
	// Alpha contains the non-randomness parameters.
	Alpha [][]float64

	// Tau defines the temperature-dependent interaction parameters.
	Tau nrtl.TauCorrelation
}

// Model constructs the corresponding NRTL activity model.
func (r NRTL) Model() activity.Model {
	return nrtl.NRTL{
		Alpha: r.Alpha,
		Tau:   r.Tau,
	}
}

// Margules supplies the parameters required to construct a
// two-component Margules activity-coefficient model.
type Margules struct {
	// A12 is the interaction parameter for component 1 wrt component 2.
	A12 float64

	// A21 is the interaction parameter for component 2 wrt component 1.
	A21 float64
}

// Model constructs the corresponding Margules activity model.
func (m Margules) Model() activity.Model {
	return margules.Margules{
		A12: m.A12,
		A21: m.A21,
	}
}

// VanLaar supplies the parameters required to construct a
// Van Laar activity-coefficient model.
type VanLaar struct {
	// A12 is the interaction parameter for component 1 wrt component 2.
	A12 float64

	// A21 is the interaction parameter for component 2 wrt component 1.
	A21 float64
}

// Model constructs the corresponding Van Laar activity model.
func (v VanLaar) Model() activity.Model {
	return vanlaar.VanLaar{
		A12: v.A12,
		A21: v.A21,
	}
}
