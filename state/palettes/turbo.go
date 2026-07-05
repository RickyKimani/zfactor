package palettes

import "image/color"

var turboColors = [...]color.Color{
	color.RGBA{R: 48, G: 18, B: 59, A: 255},   // #30123B
	color.RGBA{R: 50, G: 87, B: 220, A: 255},  // #3257DC
	color.RGBA{R: 36, G: 171, B: 220, A: 255}, // #24ABDC
	color.RGBA{R: 48, G: 214, B: 107, A: 255}, // #30D66B
	color.RGBA{R: 164, G: 252, B: 60, A: 255}, // #A4FC3C
	color.RGBA{R: 254, G: 221, B: 40, A: 255}, // #FEDD28
	color.RGBA{R: 247, G: 127, B: 12, A: 255}, // #F77F0C
	color.RGBA{R: 197, G: 37, B: 3, A: 255},   // #C52503
	color.RGBA{R: 122, G: 4, B: 3, A: 255},    // #7A0403
}

type turbo struct{}

func TurboPalette() Palette {
	return turbo{}
}

func (turbo) Isotherm(i int) color.Color {
	return turboColors[i%len(turboColors)]
}
