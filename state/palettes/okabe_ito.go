package palettes

import "github.com/rickykimani/zfactor/state/color"

var okabeItoColors = [...]color.Color{
	color.RGBA(230, 159, 0, 255),
	color.RGBA(86, 180, 233, 255),
	color.RGBA(0, 158, 115, 255),
	color.RGBA(240, 228, 66, 255),
	color.RGBA(0, 114, 178, 255),
	color.RGBA(213, 94, 0, 255),
	color.RGBA(204, 121, 167, 255),
	color.RGBA(0, 0, 0, 255),
}

type okabeIto struct{}

func OkabeItoPalette() Palette {
	return okabeIto{}
}

func (okabeIto) Isotherm(i int) color.Color {
	return okabeItoColors[i%len(okabeItoColors)]
}

func (okabeIto) CriticalIsotherm() color.Color {
	return color.Magenta
}

func (okabeIto) Dome() color.Color {
	return color.Black
}

func (okabeIto) StatePoint() color.Color {
	return color.Red
}
