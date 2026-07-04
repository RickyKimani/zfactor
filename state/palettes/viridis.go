package palettes

import "github.com/rickykimani/zfactor/state/color"

var viridisColors = [...]color.Color{
	color.RGBA(68, 1, 84, 255),    // #440154
	color.RGBA(72, 40, 120, 255),  // #482878
	color.RGBA(62, 74, 137, 255),  // #3E4A89
	color.RGBA(49, 104, 142, 255), // #31688E
	color.RGBA(38, 130, 142, 255), // #26828E
	color.RGBA(31, 158, 137, 255), // #1F9E89
	color.RGBA(53, 183, 121, 255), // #35B779
	color.RGBA(109, 205, 89, 255), // #6DCD59
	color.RGBA(180, 222, 44, 255), // #B4DE2C
	color.RGBA(253, 231, 37, 255), // #FDE725
}

type viridis struct{}

func ViridisPalette() Palette {
	return viridis{}
}

func (viridis) Isotherm(i int) color.Color {
	return viridisColors[i%len(viridisColors)]
}

func (viridis) CriticalIsotherm() color.Color {
	return color.Magenta
}

func (viridis) Dome() color.Color {
	return color.Black
}

func (viridis) StatePoint() color.Color {
	return color.Red
}
