package palettes

import "image/color"

var plasmaColors = [...]color.Color{
	color.RGBA{R: 13, G: 8, B: 135, A: 255},   // #0D0887
	color.RGBA{R: 75, G: 3, B: 161, A: 255},   // #4B03A1
	color.RGBA{R: 125, G: 3, B: 168, A: 255},  // #7D03A8
	color.RGBA{R: 168, G: 34, B: 150, A: 255}, // #A82296
	color.RGBA{R: 203, G: 70, B: 121, A: 255}, // #CB4679
	color.RGBA{R: 229, G: 107, B: 93, A: 255}, // #E56B5D
	color.RGBA{R: 248, G: 148, B: 65, A: 255}, // #F89441
	color.RGBA{R: 253, G: 195, B: 40, A: 255}, // #FDC328
	color.RGBA{R: 240, G: 249, B: 33, A: 255}, // #F0F921
}

type plasma struct{}

func PlasmaPalette() Palette {
	return plasma{}
}

func (plasma) Isotherm(i int) color.Color {
	return plasmaColors[i%len(plasmaColors)]
}
