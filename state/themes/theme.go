// Package themes provides visual themes for thermodynamic diagrams.
//
// A Theme controls the overall appearance of a plot, including its background,
// text, axes, tick marks, legend, and other user interface elements. Themes do
// not determine the colors of thermodynamic objects such as isotherms or state
// points; those are supplied by a Palette from the palettes package.
//
// Themes and palettes are designed to be orthogonal. Any palette may be used
// with any theme, allowing users to independently customize the plot's
// appearance and the colors used to represent thermodynamic data.
//
// The package provides several predefined themes.
package themes

import (
	"image/color"

	"github.com/rickykimani/zfactor/state/palettes"
	"gonum.org/v1/plot"
)

// Theme controls the appearance of a plot.
type Theme interface {
	// Apply applies the theme to the plot.
	Apply(*plot.Plot)

	// Palette returns the palette used to draw thermodynamic entities.
	Palette() palettes.Palette

	// WithPalette returns a copy of the theme using p.
	WithPalette(palettes.Palette) Theme

	// Structural/reference colors (must contrast with the theme's background)

	// Dome returns the color of the saturation dome.
	Dome() color.Color

	// CriticalIsotherm returns the color of the critical isotherm.
	CriticalIsotherm() color.Color

	// CriticalPoint returns the color of the critical point marker.
	CriticalPoint() color.Color

	// StatePoint returns the color of plotted state points.
	StatePoint() color.Color

	// IsothermLabel returns the color of the isotherm labels.
	IsothermLabel() color.Color

	// StateNumber returns the color of state point labels.
	StateNumber() color.Color

	// GridColor returns the color of grid lines.
	GridColor() color.Color
}

var (
	Default = DefaultTheme()
	Dark    = DarkTheme()
)

type defaultTheme struct {
	palette palettes.Palette
}

func (defaultTheme) Apply(p *plot.Plot) {
	// Title
	p.Title.TextStyle.Color = color.Black

	// Axis labels
	p.X.Label.TextStyle.Color = color.Black
	p.Y.Label.TextStyle.Color = color.Black

	// Axis lines
	p.X.LineStyle.Color = color.Black
	p.Y.LineStyle.Color = color.Black

	// Tick marks
	p.X.Tick.LineStyle.Color = color.Black
	p.Y.Tick.LineStyle.Color = color.Black

	// Tick labels
	p.X.Tick.Label.Color = color.Black
	p.Y.Tick.Label.Color = color.Black

	// Legend
	p.Legend.TextStyle.Color = color.Black
}

func (t defaultTheme) Palette() palettes.Palette {
	return t.palette
}

func (t defaultTheme) WithPalette(p palettes.Palette) Theme {
	t.palette = p
	return t
}

func (defaultTheme) Dome() color.Color {
	return color.Black
}

func (defaultTheme) CriticalIsotherm() color.Color {
	return color.RGBA{R: 255, G: 0, B: 255, A: 255} // Magenta
}

func (defaultTheme) CriticalPoint() color.Color {
	return color.RGBA{R: 255, G: 0, B: 255, A: 255} // Magenta
}

func (defaultTheme) StatePoint() color.Color {
	return color.RGBA{R: 255, G: 0, B: 0, A: 255} // Red
}

func (defaultTheme) IsothermLabel() color.Color {
	return color.Black
}

func (defaultTheme) StateNumber() color.Color {
	return color.Black
}

func (defaultTheme) GridColor() color.Color {
	return color.RGBA{R: 220, G: 220, B: 220, A: 255} // Light grey
}

// DefaultTheme provides the standard light theme as Gonum defaults to a white background
func DefaultTheme() Theme {
	return defaultTheme{
		palette: palettes.DefaultPalette(),
	}
}
