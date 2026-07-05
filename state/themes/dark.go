package themes

import (
	"image/color"

	"github.com/rickykimani/zfactor/state/palettes"
	"gonum.org/v1/plot"
)

type darkTheme struct {
	palette palettes.Palette
}

func DarkTheme() Theme {
	return darkTheme{palette: palettes.DefaultPalette()}
}

func (t darkTheme) Apply(p *plot.Plot) {
	// Background
	p.BackgroundColor = color.RGBA{R: 24, G: 24, B: 24, A: 255}

	foreground := color.RGBA{R: 230, G: 230, B: 230, A: 255}

	// Title
	p.Title.TextStyle.Color = foreground

	// Axis labels
	p.X.Label.TextStyle.Color = foreground
	p.Y.Label.TextStyle.Color = foreground

	// Axis lines
	p.X.LineStyle.Color = foreground
	p.Y.LineStyle.Color = foreground

	// Tick marks
	p.X.Tick.LineStyle.Color = foreground
	p.Y.Tick.LineStyle.Color = foreground

	// Tick labels
	p.X.Tick.Label.Color = foreground
	p.Y.Tick.Label.Color = foreground

	// Legend
	p.Legend.TextStyle.Color = foreground
}

func (t darkTheme) Palette() palettes.Palette {
	return t.palette
}

func (t darkTheme) WithPalette(p palettes.Palette) Theme {
	t.palette = p
	return t
}

func (darkTheme) Dome() color.Color {
	return color.RGBA{R: 200, G: 200, B: 200, A: 255} // Light grey
}

func (darkTheme) CriticalIsotherm() color.Color {
	return color.RGBA{R: 255, G: 100, B: 255, A: 255} // Bright magenta
}

func (darkTheme) CriticalPoint() color.Color {
	return color.RGBA{R: 255, G: 100, B: 255, A: 255} // Bright magenta
}

func (darkTheme) StatePoint() color.Color {
	return color.RGBA{R: 255, G: 80, B: 80, A: 255} // Bright light red
}

func (darkTheme) IsothermLabel() color.Color {
	return color.RGBA{R: 180, G: 180, B: 180, A: 255} // Medium grey
}

func (darkTheme) StateNumber() color.Color {
	return color.RGBA{R: 180, G: 180, B: 180, A: 255} // Medium grey
}

func (darkTheme) GridColor() color.Color {
	return color.RGBA{R: 50, G: 50, B: 50, A: 255} // Dark grey
}
