package themes

import (
	"github.com/rickykimani/zfactor/state/color"
	"github.com/rickykimani/zfactor/state/palettes"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/vg"
)

type darkTheme struct {
	palette palettes.Palette
}

func DarkTheme() Theme {
	return darkTheme{palette: palettes.DefaultPalette()}
}

func (t darkTheme) Apply(p *plot.Plot) {
	// Background
	p.BackgroundColor = color.RGBA(24, 24, 24, 255)

	foreground := color.RGBA(230, 230, 230, 255)

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
	p.Legend.Top = true
	p.Legend.Left = false

	// Optional border around legend.
	p.Legend.Padding = vg.Points(4)
}

func (t darkTheme) Palette() palettes.Palette {
	return t.palette
}

func (t darkTheme) IsothermLabel() color.Color {
	return color.Grey
}

func (t darkTheme) StateNumber() color.Color {
	return color.Grey
}

func (t darkTheme) WithPalette(p palettes.Palette) Theme {
	t.palette = p
	return t
}
