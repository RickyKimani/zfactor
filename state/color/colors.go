// Package color provides reusable color definitions and utilities for
// thermodynamic diagram rendering.
//
// The package defines a Color type which is an alias of Go's standard image/color
// interfaces together with a collection of commonly used named colors.
//
// It is primarily intended for use by the palettes and themes packages, but
// may also be used directly when defining custom plot styles.
package color

import "image/color"

// Color is an alias for image/color.Color, representing colors in the plot.
type Color = color.Color

func RGBA(r, g, b, a uint8) Color {
	return color.RGBA{R: r, G: g, B: b, A: a}
}

// Standard colors provided for convenience.
var (
	Red     Color = color.RGBA{R: 255, G: 0, B: 0, A: 255}
	Green   Color = color.RGBA{R: 0, G: 255, B: 0, A: 255}
	Blue    Color = color.RGBA{R: 0, G: 0, B: 255, A: 255}
	Yellow  Color = color.RGBA{R: 255, G: 255, B: 0, A: 255}
	Cyan    Color = color.RGBA{R: 0, G: 255, B: 255, A: 255}
	Magenta Color = color.RGBA{R: 255, G: 0, B: 255, A: 255}
	White   Color = color.RGBA{R: 255, G: 255, B: 255, A: 255}
	Black   Color = color.RGBA{R: 0, G: 0, B: 0, A: 255}
	Pink    Color = color.RGBA{R: 255, G: 192, B: 203, A: 255}
	Orange  Color = color.RGBA{R: 255, G: 165, B: 0, A: 255}
	Purple  Color = color.RGBA{R: 128, G: 0, B: 128, A: 255}
	Grey    Color = color.RGBA{R: 128, G: 128, B: 128, A: 255}
)
