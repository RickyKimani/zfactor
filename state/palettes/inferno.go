package palettes

import "image/color"

var infernoColors = [...]color.Color{
	color.RGBA{R: 0, G: 0, B: 4, A: 255},       // #000004
	color.RGBA{R: 31, G: 12, B: 72, A: 255},    // #1F0C48
	color.RGBA{R: 85, G: 15, B: 109, A: 255},   // #550F6D
	color.RGBA{R: 136, G: 34, B: 106, A: 255},  // #88226A
	color.RGBA{R: 186, G: 54, B: 85, A: 255},   // #BA3655
	color.RGBA{R: 227, G: 89, B: 51, A: 255},   // #E35933
	color.RGBA{R: 249, G: 140, B: 10, A: 255},  // #F98C0A
	color.RGBA{R: 249, G: 201, B: 50, A: 255},  // #F9C932
	color.RGBA{R: 252, G: 255, B: 164, A: 255}, // #FCFFA4
}

type inferno struct{}

func InfernoPalette() Palette {
	return inferno{}
}

func (inferno) Isotherm(i int) color.Color {
	return infernoColors[i%len(infernoColors)]
}
