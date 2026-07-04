package palettes

import "github.com/rickykimani/zfactor/state/color"

var turboColors = [...]color.Color{
	color.RGBA(48, 18, 59, 255),
	color.RGBA(50, 87, 220, 255),
	color.RGBA(36, 171, 220, 255),
	color.RGBA(48, 214, 107, 255),
	color.RGBA(164, 252, 60, 255),
	color.RGBA(254, 221, 40, 255),
	color.RGBA(247, 127, 12, 255),
	color.RGBA(197, 37, 3, 255),
	color.RGBA(122, 4, 3, 255),
}

type turbo struct{}

func TurboPalette() Palette {
	return turbo{}
}

func (turbo) Isotherm(i int) color.Color {
	return turboColors[i%len(turboColors)]
}

func (turbo) CriticalIsotherm() color.Color {
	return color.Magenta
}

func (turbo) Dome() color.Color {
	return color.Black
}

func (turbo) StatePoint() color.Color {
	return color.Red
}
