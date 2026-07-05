package palettes

import "image/color"

var viridisColors = [...]color.Color{
	color.RGBA{R: 68, G: 1, B: 84, A: 255},    // #440154
	color.RGBA{R: 72, G: 40, B: 120, A: 255},  // #482878
	color.RGBA{R: 62, G: 74, B: 137, A: 255},  // #3E4A89
	color.RGBA{R: 49, G: 104, B: 142, A: 255}, // #31688E
	color.RGBA{R: 38, G: 130, B: 142, A: 255}, // #26828E
	color.RGBA{R: 31, G: 158, B: 137, A: 255}, // #1F9E89
	color.RGBA{R: 53, G: 183, B: 121, A: 255}, // #35B779
	color.RGBA{R: 109, G: 205, B: 89, A: 255}, // #6DCD59
	color.RGBA{R: 180, G: 222, B: 44, A: 255}, // #B4DE2C
	color.RGBA{R: 253, G: 231, B: 37, A: 255}, // #FDE725
}

type viridis struct{}

func ViridisPalette() Palette {
	return viridis{}
}

func (viridis) Isotherm(i int) color.Color {
	return viridisColors[i%len(viridisColors)]
}
