package palettes

import "image/color"

var okabeItoColors = [...]color.Color{
	color.RGBA{R: 230, G: 159, B: 0, A: 255},   // #E69F00
	color.RGBA{R: 86, G: 180, B: 233, A: 255},  // #56B4E9
	color.RGBA{R: 0, G: 158, B: 115, A: 255},   // #009E73
	color.RGBA{R: 240, G: 228, B: 66, A: 255},  // #F0E442
	color.RGBA{R: 0, G: 114, B: 178, A: 255},   // #0072B2
	color.RGBA{R: 213, G: 94, B: 0, A: 255},    // #D55E00
	color.RGBA{R: 204, G: 121, B: 167, A: 255}, // #CC79A7
	color.RGBA{R: 0, G: 0, B: 0, A: 255},       // #000000
}

type okabeIto struct{}

func OkabeItoPalette() Palette {
	return okabeIto{}
}

func (okabeIto) Isotherm(i int) color.Color {
	return okabeItoColors[i%len(okabeItoColors)]
}
