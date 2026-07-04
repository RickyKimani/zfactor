package palettes

import "github.com/rickykimani/zfactor/state/color"

var plasmaColors = [...]color.Color{
	color.RGBA(13, 8, 135, 255),
	color.RGBA(75, 3, 161, 255),
	color.RGBA(125, 3, 168, 255),
	color.RGBA(168, 34, 150, 255),
	color.RGBA(203, 70, 121, 255),
	color.RGBA(229, 107, 93, 255),
	color.RGBA(248, 148, 65, 255),
	color.RGBA(253, 195, 40, 255),
	color.RGBA(240, 249, 33, 255),
}

type plasma struct{}

func PlasmaPalette() Palette {
	return plasma{}
}

func (plasma) Isotherm(i int) color.Color {
	return plasmaColors[i%len(plasmaColors)]
}

func (plasma) CriticalIsotherm() color.Color {
	return color.Magenta
}

func (plasma) Dome() color.Color {
	return color.Black
}

func (plasma) StatePoint() color.Color {
	return color.Red
}
