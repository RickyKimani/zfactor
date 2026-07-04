package state

import "gonum.org/v1/plot/vg"

// Length is an alias for vg.Length, representing physical length units for plotting.
type Length = vg.Length

// Common length units for specifying plot dimensions.
const (
	Inch       Length = vg.Inch
	Centimeter Length = vg.Centimeter
	Millimeter Length = vg.Millimeter
)
