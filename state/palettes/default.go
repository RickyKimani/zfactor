package palettes

import (
	"github.com/rickykimani/zfactor/state/color"
)

type defaultPalette struct{}

func DefaultPalette() Palette {
	return defaultPalette{}
}

func (defaultPalette) Isotherm(int) color.Color {
	return color.Blue
}

func (defaultPalette) CriticalIsotherm() color.Color {
	return color.Magenta
}

func (defaultPalette) Dome() color.Color {
	return color.Black
}

func (defaultPalette) StatePoint() color.Color {
	return color.Red
}
