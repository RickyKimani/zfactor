package palettes

import (
	"image/color"
	"github.com/rickykimani/zfactor/state/statecolor"
)

type defaultPalette struct{}

func DefaultPalette() Palette {
	return defaultPalette{}
}

func (defaultPalette) Isotherm(int) color.Color {
	return statecolor.Blue // #0000FF
}
