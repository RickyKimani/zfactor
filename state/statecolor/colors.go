// Package statecolor provides standard color variables using the standard library image/color types.
package statecolor

import "image/color"

// Standard colors provided for convenience.
var (
	Red     = color.RGBA{R: 255, G: 0, B: 0, A: 255}
	Green   = color.RGBA{R: 0, G: 255, B: 0, A: 255}
	Blue    = color.RGBA{R: 0, G: 0, B: 255, A: 255}
	Yellow  = color.RGBA{R: 255, G: 255, B: 0, A: 255}
	Cyan    = color.RGBA{R: 0, G: 255, B: 255, A: 255}
	Magenta = color.RGBA{R: 255, G: 0, B: 255, A: 255}
	White   = color.RGBA{R: 255, G: 255, B: 255, A: 255}
	Black   = color.RGBA{R: 0, G: 0, B: 0, A: 255}
	Pink    = color.RGBA{R: 255, G: 192, B: 203, A: 255}
	Orange  = color.RGBA{R: 255, G: 165, B: 0, A: 255}
	Purple  = color.RGBA{R: 128, G: 0, B: 128, A: 255}
	Grey    = color.RGBA{R: 128, G: 128, B: 128, A: 255}
)
