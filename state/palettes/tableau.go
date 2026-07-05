package palettes

import "image/color"

var tableauColors = [...]color.Color{
	color.RGBA{R: 78, G: 121, B: 167, A: 255},  // #4E79A7
	color.RGBA{R: 242, G: 142, B: 43, A: 255},  // #F28E2B
	color.RGBA{R: 225, G: 87, B: 89, A: 255},   // #E15759
	color.RGBA{R: 118, G: 183, B: 178, A: 255}, // #76B7B2
	color.RGBA{R: 89, G: 161, B: 79, A: 255},   // #59A14F
	color.RGBA{R: 237, G: 201, B: 72, A: 255},  // #EDC948
	color.RGBA{R: 176, G: 122, B: 161, A: 255}, // #B07AA1
	color.RGBA{R: 255, G: 157, B: 167, A: 255}, // #FF9DA7
	color.RGBA{R: 156, G: 117, B: 95, A: 255},  // #9C755F
	color.RGBA{R: 186, G: 176, B: 172, A: 255}, // #BAB0AC
}

type tableau struct{}

func TableauPalette() Palette {
	return tableau{}
}

func (tableau) Isotherm(i int) color.Color {
	return tableauColors[i%len(tableauColors)]
}
