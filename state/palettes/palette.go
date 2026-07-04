// Package palettes provides reusable color palettes for thermodynamic
// diagrams.
//
// A Palette defines the colors used to render thermodynamic entities such as
// isotherms, saturation domes, critical isotherms, and state points.
//
// Palettes are independent of the overall appearance of the plot. They are
// intended to be used together with a theme from the themes package, allowing
// the same palette to be displayed using different visual styles (for example,
// light or dark backgrounds).
//
// The package includes several built-in palettes suitable for different use cases.
package palettes

import "github.com/rickykimani/zfactor/state/color"

var (
	Default  = DefaultPalette()
	Viridis  = ViridisPalette()
	OkabeIto = OkabeItoPalette()
	Turbo    = TurboPalette()
	Inferno  = InfernoPalette()
	Plasma   = PlasmaPalette()
	Tableau  = TableauPalette()
)

// Palette defines the colors used to render thermodynamic entities
// on a state diagram.
type Palette interface {
	// Isotherm returns the color of the i-th isotherm.
	Isotherm(i int) color.Color

	// CriticalIsotherm returns the color of the critical isotherm.
	CriticalIsotherm() color.Color

	// Dome returns the color of the saturation dome.
	Dome() color.Color

	// StatePoint returns the color of plotted state points.
	StatePoint() color.Color
}
