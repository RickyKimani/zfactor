// Package state provides functionality for defining thermodynamic states and generating
// visual representations such as PV diagrams.
package state

import (
	"errors"
	"fmt"
	"image/color"
	"os"
	"path/filepath"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/substance"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

var validExts = map[string]bool{
	".eps":  true,
	".jpg":  true,
	".jpeg": true,
	".pdf":  true,
	".png":  true,
	".svg":  true,
	".tex":  true,
	".tif":  true,
	".tiff": true,
}

// Color is an alias for image/color.Color, representing colors in the plot.
type Color = color.Color

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

// Length is an alias for vg.Length, representing physical length units for plotting.
type Length = vg.Length

// Common length units for specifying plot dimensions.
const (
	Inch       Length = vg.Inch
	Centimeter Length = vg.Centimeter
	Millimeter Length = vg.Millimeter
)

// State represents a specific thermodynamic state of a substance defined by its
// temperature and pressure.
type State struct {
	Substance   *substance.Substance
	Temperature float64 // Temperature in Kelvin
	Pressure    float64 // Pressure in bar
}

// NewState creates a new State object. It validates that the temperature and pressure
// are positive values.
func NewState(substance *substance.Substance, t, p float64) (*State, error) {
	if t <= 0 {
		return nil, zfactor.ErrTemp
	}
	if p <= 0 {
		return nil, zfactor.ErrPressure
	}
	return &State{
		Substance:   substance,
		Temperature: t,
		Pressure:    p,
	}, nil
}

// PVConfig holds configuration options for customizing the appearance of the PV diagram.
type PVConfig struct {
	// Type specifies the cubic Equation of State (EOS) model to use for generating the PV diagram.
	// This field is required; DrawPV will return an error if it is nil.
	Type cubic.EOSType
	// Title is the title of the plot. If empty, a default title is generated.
	Title string
	// TitleColor is the color of the title text. Defaults to black if nil.
	TitleColor Color
	// XLabelColor is the color of the X axis label text. Defaults to black if nil
	XLabelColor Color
	// YLabelColor is the color of the Y axis label text. Defaults to black if nil
	YLabelColor Color
	// Width is the width of the output image. Defaults to 6 inches if 0.
	Width Length
	// Height is the height of the output image. Defaults to 4 inches if 0.
	Height Length
	// IsothermsColor is the color of the isotherm lines. Defaults to blue if nil.
	IsothermsColor Color
	// CriticalIsothermColor is the color of the critical isotherm (T=Tc). Defaults to magenta if nil.
	CriticalIsothermColor Color
	// DomeColor is the color of the saturation dome. Defaults to black if nil.
	DomeColor Color
	// StatePointColor is the color of the point representing the state. Defaults to red if nil.
	StatePointColor Color
	// NumberStates places a number on alongside the state point in the order they occur in states ...*State
	NumberStates bool
	// StatePointNumberColor is the color of the number of the state. Defaults to black if nil.
	StatePointNumberColor Color
	// LabelIsotherms places a label alongside the isotherm with the numerical value of the temperature
	LabelIsotherms bool
	// IsothermLabelColor is the color of the isotherm label. Defaults to black if nil.
	IsothermLabelColor Color
	// VolumeScaleFactor determines the maximum volume shown on the X-axis as a multiple of the critical volume (Vc).
	// If 0, it defaults to 7.0.
	VolumeScaleFactor float64
	// ShowOutputPath determines whether to print the full path of the saved image to stdout upon success.
	ShowOutputPath bool
}

// DrawPV generates a Pressure-Volume (PV) diagram for the provided states.
// It plots the critical isotherm, the saturation dome (two-phase region), and the
// specific isotherms for each state provided. The resulting plot is saved to the
// file specified by 'output'.
func DrawPV(cfg *PVConfig, output string, states ...*State) error {
	if cfg == nil {
		return errors.New("configuration error: config cannot be nil")
	}
	if cfg.Type == nil {
		return errors.New("configuration error: 'Type' field (EOS model) is required")
	}
	ext := filepath.Ext(output)
	if ok := validExts[ext]; !ok {
		closest := ""
		minDist := int(^uint(0) >> 1)
		for valid := range validExts {
			dist := levenshtein(ext, valid)
			if dist < minDist {
				minDist = dist
				closest = valid
			}
		}
		suggestion := output[:len(output)-len(ext)] + closest
		return fmt.Errorf("invalid file extension: %s. Did you mean %q instead?", output, suggestion)
	}
	name, err := verifySubstances(states...)
	if err != nil {
		return fmt.Errorf("oops, something went wrong: %w", err)
	}
	p := plot.New()

	if cfg.Title == "" {
		p.Title.Text = fmt.Sprintf("PV Diagram for %s", name)
	} else {
		p.Title.Text = cfg.Title
	}

	if cfg.TitleColor != nil {
		p.Title.TextStyle.Color = cfg.TitleColor
	}

	p.X.Label.Text = "Molar Volume (cm³/mol)"
	if cfg.XLabelColor != nil {
		p.X.Label.TextStyle.Color = cfg.XLabelColor
	}
	p.Y.Label.Text = "Pressure (bar)"
	if cfg.YLabelColor != nil {
		p.X.Label.TextStyle.Color = cfg.YLabelColor
	}

	// Use Linear Scale but be smart about limits
	// p.X.Scale = plot.LogScale{}

	const R = zfactor.RSI * 10 // bar*cm^3/(mol*K)

	s0 := states[0]
	Tc := s0.Substance.Critical.Tc
	Pc := s0.Substance.Critical.Pc
	Vc := s0.Substance.Critical.Vc

	// 1. Draw Critical Isotherm (T = Tc)
	// This defines the boundary between subcritical and supercritical
	critCfg := s0.Substance.CubicConfig(cfg.Type, Tc, Pc, R)
	b := critCfg.Type.Params().Omega * R * Tc / Pc

	// Define V range based on Vc
	// Start near b, go up to a reasonable multiple of Vc
	minV := b * 1.1
	// Default max view: if Vc is known, use it. Else guess.
	maxViewV := minV * 15
	if Vc > 0 {
		factor := cfg.VolumeScaleFactor
		if factor <= 0 {
			factor = 7.0
		}
		maxViewV = Vc * factor
	}

	// Check if any state is outside this view
	for _, s := range states {
		// Estimate V for state
		estV := R * s.Temperature / s.Pressure
		if estV > maxViewV {
			maxViewV = estV * 1.1
		}
	}

	critPts := make(plotter.XYs, 0)
	// Generate points for Critical Isotherm
	// Use logarithmic spacing for smoothness even on linear plot
	for v := minV; v <= maxViewV; v *= 1.05 {
		presRes, err := cubic.Pressure(critCfg, v)
		if err == nil && presRes.P > 0 {
			critPts = append(critPts, plotter.XY{X: v, Y: presRes.P})
		}
	}
	critLine, _ := plotter.NewLine(critPts)
	if cfg.CriticalIsothermColor == nil {
		critLine.Color = Magenta
	} else {
		critLine.Color = cfg.CriticalIsothermColor
	}
	critLine.LineStyle.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
	critLine.LineStyle.Width = vg.Points(1)
	p.Add(critLine)

	if cfg.LabelIsotherms && len(critPts) > 0 {
		lastPt := critPts[len(critPts)-1]
		labels, _ := plotter.NewLabels(plotter.XYLabels{
			XYs:    []plotter.XY{lastPt},
			Labels: []string{fmt.Sprintf("Tc=%.1f K", Tc)},
		})
		labels.Offset.X = vg.Points(2)
		if cfg.IsothermLabelColor != nil {
			labels.TextStyle[0].Color = cfg.IsothermLabelColor
		}
		p.Add(labels)
	}

	// 2. Draw Saturation Dome
	domeCfg := s0.Substance.CubicConfig(cfg.Type, Tc, Pc, R)
	var liquidPts, vaporPts plotter.XYs

	// Range from 0.6 Tc to 0.99 Tc
	// Closer to Tc is harder to converge
	startT := Tc * 0.6
	endT := Tc * 0.99
	stepT := (endT - startT) / 100

	for t := startT; t <= endT; t += stepT {
		pSat, err := cubic.SaturationPressure(domeCfg, t)
		if err != nil {
			continue
		}
		domeCfg.T = t
		domeCfg.P = pSat
		volRes, err := cubic.SolveForVolume(domeCfg)
		if err != nil {
			continue
		}
		roots := volRes.Clean()
		if len(roots) >= 2 {
			liquidPts = append(liquidPts, plotter.XY{X: roots[0], Y: pSat})
			vaporPts = append(vaporPts, plotter.XY{X: roots[len(roots)-1], Y: pSat})
		}
	}

	// Add Critical Point to close the dome
	if Vc > 0 {
		liquidPts = append(liquidPts, plotter.XY{X: Vc, Y: Pc})
	}

	// Connect vapor points back to liquid (reverse order)
	for i := len(vaporPts) - 1; i >= 0; i-- {
		liquidPts = append(liquidPts, vaporPts[i])
	}

	if len(liquidPts) > 0 {
		domeLine, _ := plotter.NewLine(liquidPts)
		if cfg.DomeColor == nil {
			domeLine.Color = Black
		} else {
			domeLine.Color = cfg.DomeColor
		}
		domeLine.LineStyle.Width = vg.Points(1.5)
		p.Add(domeLine)
	}

	// 3. Mark Critical Point
	if Vc > 0 {
		cp, _ := plotter.NewScatter(plotter.XYs{{X: Vc, Y: Pc}})
		cp.GlyphStyle.Shape = draw.CrossGlyph{}
		cp.Color = color.RGBA{R: 0, A: 255}
		p.Add(cp)
	}

	// 4. Draw States and their Isotherms
	for i, state := range states {
		stateCfg := state.Substance.CubicConfig(cfg.Type, state.Temperature, state.Pressure, R)

		// Draw Isotherm
		isoPts := make(plotter.XYs, 0)
		for v := minV; v <= maxViewV; v *= 1.05 {
			presRes, err := cubic.Pressure(stateCfg, v)
			if err == nil && presRes.P > 0 {
				isoPts = append(isoPts, plotter.XY{X: v, Y: presRes.P})
			}
		}
		isoLine, _ := plotter.NewLine(isoPts)
		if cfg.IsothermsColor == nil {
			isoLine.Color = Blue
		} else {
			isoLine.Color = cfg.IsothermsColor
		}
		p.Add(isoLine)

		if cfg.LabelIsotherms && len(isoPts) > 0 {
			lastPt := isoPts[len(isoPts)-1]
			labels, _ := plotter.NewLabels(plotter.XYLabels{
				XYs:    []plotter.XY{lastPt},
				Labels: []string{fmt.Sprintf("T=%.1f K", state.Temperature)},
			})
			labels.Offset.X = vg.Points(2)
			// Shift label to avoid overlap with Critical Isotherm
			if state.Temperature < Tc {
				labels.Offset.Y = vg.Points(-10)
			} else {
				labels.Offset.Y = vg.Points(10)
			}
			if cfg.IsothermLabelColor != nil {
				labels.TextStyle[0].Color = cfg.IsothermLabelColor
			}
			p.Add(labels)
		}

		// Calculate State Point
		volRes, err := cubic.SolveForVolume(stateCfg)
		if err != nil {
			continue
		}
		roots := volRes.Clean()

		if len(roots) == 0 {
			continue
		}

		// Determine which root represents the state
		// If 1 root: Supercritical or Single Phase
		// If 3 roots: Two-Phase region possible.
		// But we are given P and T.
		// If T < Tc and P < Psat -> Vapor (largest root)
		// If T < Tc and P > Psat -> Liquid (smallest root)
		// If T > Tc -> Single root

		var stateV float64

		if state.Temperature >= Tc {
			stateV = roots[0] // Only 1 real root usually
		} else {
			// Subcritical
			pSat, err := cubic.SaturationPressure(stateCfg, state.Temperature)
			if err == nil {
				if state.Pressure > pSat {
					stateV = roots[0] // Liquid
				} else if state.Pressure < pSat {
					stateV = roots[len(roots)-1] // Vapor
				} else {
					// Saturation
					// Ambiguous V, could be anywhere.
					// Usually user implies one, but let's pick Vapor for visualization or both?
					stateV = roots[len(roots)-1]
				}
			} else {
				// Fallback
				stateV = roots[len(roots)-1]
			}
		}

		// Plot State Marker
		scatter, _ := plotter.NewScatter(plotter.XYs{{X: stateV, Y: state.Pressure}})
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.GlyphStyle.Radius = vg.Points(4)
		if cfg.StatePointColor == nil {
			scatter.Color = Red
		} else {
			scatter.Color = cfg.StatePointColor
		}
		p.Add(scatter)

		if cfg.NumberStates {
			labels, _ := plotter.NewLabels(plotter.XYLabels{
				XYs:    []plotter.XY{{X: stateV, Y: state.Pressure}},
				Labels: []string{fmt.Sprintf("%d", i+1)},
			})
			labels.Offset.X = vg.Points(5)
			labels.Offset.Y = vg.Points(5)
			if cfg.StatePointNumberColor != nil {
				labels.TextStyle[0].Color = cfg.StatePointNumberColor
			}
			p.Add(labels)
		}
	}

	// Set Axes Limits
	p.X.Min = 0
	p.X.Max = maxViewV
	p.Y.Min = 0
	p.Y.Max = Pc * 1.5
	if states[0].Pressure > p.Y.Max {
		p.Y.Max = states[0].Pressure * 1.1
	}

	width := cfg.Width
	if width == 0 {
		width = 6 * vg.Inch
	}
	height := cfg.Height
	if height == 0 {
		height = 4 * vg.Inch
	}

	err = p.Save(width, height, output)
	if err != nil {
		return err
	}

	if cfg.ShowOutputPath {
		wd, err := os.Getwd()
		if err != nil {
			return fmt.Errorf("failed to get working directory: %w", err)
		}
		fmt.Printf("image saved to %s\n", filepath.Join(wd, output))
	}

	return nil
}

// verifySubstances ensures that all provided states belong to the same substance.
// It returns the name of the substance if consistent, or an error otherwise.
func verifySubstances(states ...*State) (string, error) {
	var prev string
	var curr string
	prev = states[0].Substance.Name
	for _, state := range states {
		curr = state.Substance.Name
		if curr != prev {
			return "", errors.New("substance mismatch")
		}
		prev = curr
	}
	return curr, nil
}

func levenshtein(s1, s2 string) int {
	r1, r2 := []rune(s1), []rune(s2)
	n, m := len(r1), len(r2)
	if n == 0 {
		return m
	}
	if m == 0 {
		return n
	}
	row := make([]int, n+1)
	for i := 0; i <= n; i++ {
		row[i] = i
	}
	for j := 1; j <= m; j++ {
		prev := j
		for i := 1; i <= n; i++ {
			cost := 0
			if r1[i-1] != r2[j-1] {
				cost = 1
			}
			current := min(row[i]+1, prev+1, row[i-1]+cost)
			row[i-1] = prev
			prev = current
		}
		row[n] = prev
	}
	return row[n]
}
