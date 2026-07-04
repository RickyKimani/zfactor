package palettes

import "github.com/rickykimani/zfactor/state/color"

var infernoColors = [...]color.Color{
	color.RGBA(0, 0, 4, 255),
	color.RGBA(31, 12, 72, 255),
	color.RGBA(85, 15, 109, 255),
	color.RGBA(136, 34, 106, 255),
	color.RGBA(186, 54, 85, 255),
	color.RGBA(227, 89, 51, 255),
	color.RGBA(249, 140, 10, 255),
	color.RGBA(249, 201, 50, 255),
	color.RGBA(252, 255, 164, 255),
}

type inferno struct{}

func InfernoPalette() Palette {
	return inferno{}
}

func (inferno) Isotherm(i int) color.Color {
	return infernoColors[i%len(infernoColors)]
}

func (inferno) CriticalIsotherm() color.Color {
	return color.Magenta
}

func (inferno) Dome() color.Color {
	return color.Black
}

func (inferno) StatePoint() color.Color {
	return color.Red
}
