package palettes

import "github.com/rickykimani/zfactor/state/color"

var tableauColors = [...]color.Color{
	color.RGBA(78, 121, 167, 255),
	color.RGBA(242, 142, 43, 255),
	color.RGBA(225, 87, 89, 255),
	color.RGBA(118, 183, 178, 255),
	color.RGBA(89, 161, 79, 255),
	color.RGBA(237, 201, 72, 255),
	color.RGBA(176, 122, 161, 255),
	color.RGBA(255, 157, 167, 255),
	color.RGBA(156, 117, 95, 255),
	color.RGBA(186, 176, 172, 255),
}

type tableau struct{}

func TableauPalette() Palette {
	return tableau{}
}

func (tableau) Isotherm(i int) color.Color {
	return tableauColors[i%len(tableauColors)]
}

func (tableau) CriticalIsotherm() color.Color {
	return color.Magenta
}

func (tableau) Dome() color.Color {
	return color.Black
}

func (tableau) StatePoint() color.Color {
	return color.Red
}
