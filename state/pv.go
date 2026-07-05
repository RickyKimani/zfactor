package state

import (
	"errors"
	"fmt"
	"os"
	"path/filepath"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/state/themes"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)


// PVConfig holds configuration options for customizing the appearance of the PV diagram.
type PVConfig struct {
	// Type specifies the cubic Equation of State (EOS) model to use for generating the PV diagram.
	// This field is required; DrawPV will return an error if it is nil.
	Type cubic.EOSType

	// Title is the title of the plot. If empty, a default title is generated.
	Title string

	// Theme is the theme of the plot. If unset, the default theme is used.
	Theme themes.Theme

	// Width is the width of the output image. Defaults to 6 inches if 0.
	Width Length
	// Height is the height of the output image. Defaults to 4 inches if 0.
	Height Length

	// NumberStates places a number on alongside the state point in the order they occur in states ...*State
	NumberStates bool

	// LabelIsotherms places a label alongside the isotherm with the numerical value of the temperature
	LabelIsotherms bool

	// VolumeScaleFactor determines the maximum volume shown on the X-axis as a multiple of the critical volume (Vc).
	// If 0, it defaults to 7.0.
	VolumeScaleFactor float64

	// ShowOutputPath determines whether to print the full path of the saved image to stdout upon success.
	ShowOutputPath bool
}

func DefaultPVConfig(eos cubic.EOSType) *PVConfig {
	return &PVConfig{
		Type:              eos,
		Theme:             themes.DefaultTheme(),
		Width:             6 * Inch,
		Height:            4 * Inch,
		VolumeScaleFactor: 7,
		ShowOutputPath:    true,
	}
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
	theme := cfg.Theme
	if cfg.Theme == nil {
		theme = themes.DefaultTheme()
	}
	palette := theme.Palette()
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

	theme.Apply(p)

	// Add grid lines
	grid := plotter.NewGrid()
	grid.Horizontal.Color = theme.GridColor()
	grid.Vertical.Color = theme.GridColor()
	p.Add(grid)

	if cfg.Title == "" {
		p.Title.Text = fmt.Sprintf("PV Diagram for %s", name)
	} else {
		p.Title.Text = cfg.Title
	}

	p.X.Label.Text = "Molar Volume (cm³/mol)"
	p.Y.Label.Text = "Pressure (bar)"

	// Use Linear Scale but be smart about limits
	// p.X.Scale = plot.LogScale{}

	const R = zfactor.RSI * 10 // bar*cm^3/(mol*K)

	s0 := states[0]
	Tc := s0.Substance.Critical.Tc
	Pc := s0.Substance.Critical.Pc
	Vc := s0.Substance.Critical.Vc

	// 1. Draw Critical Isotherm (T = Tc)
	// This defines the boundary between subcritical and supercritical
	critCfg := s0.Substance.CubicConfig(cfg.Type, zfactor.Args{T: Tc, P: Pc, R: R})
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
	critLine, err := plotter.NewLine(critPts)
	if err != nil {
		return fmt.Errorf("oops, something went wrong: %w", err)
	}
	critLine.Color = theme.CriticalIsotherm()

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
		labels.TextStyle[0].Color = theme.IsothermLabel()
		p.Add(labels)
	}

	// 2. Draw Saturation Dome
	domeCfg := s0.Substance.CubicConfig(cfg.Type, zfactor.Args{T: Tc, P: Pc, R: R})
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
		domeLine, err := plotter.NewLine(liquidPts)
		if err != nil {
			return fmt.Errorf("oops, something went wrong: %w", err)
		}
		domeLine.Color = theme.Dome()
		domeLine.LineStyle.Width = vg.Points(1.5)
		p.Add(domeLine)
	}

	// 3. Mark Critical Point
	if Vc > 0 {
		cp, _ := plotter.NewScatter(plotter.XYs{{X: Vc, Y: Pc}})
		cp.GlyphStyle.Shape = draw.CrossGlyph{}
		cp.Color = theme.CriticalPoint()
		p.Add(cp)
	}

	// 4. Draw States and their Isotherms
	for i, state := range states {
		stateCfg := state.Substance.CubicConfig(cfg.Type, zfactor.Args{T: state.Temperature, P: state.Pressure, R: R})

		// Draw Isotherm
		isoPts := make(plotter.XYs, 0)
		for v := minV; v <= maxViewV; v *= 1.05 {
			presRes, err := cubic.Pressure(stateCfg, v)
			if err == nil && presRes.P > 0 {
				isoPts = append(isoPts, plotter.XY{X: v, Y: presRes.P})
			}
		}
		isoLine, err := plotter.NewLine(isoPts)
		if err != nil {
			return fmt.Errorf("oops, something went wrong: %w", err)
		}
		isoLine.Color = palette.Isotherm(i)
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
			labels.TextStyle[0].Color = theme.IsothermLabel()

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
		scatter, err := plotter.NewScatter(plotter.XYs{{X: stateV, Y: state.Pressure}})
		if err != nil {
			return fmt.Errorf("oops, something went wrong: %w", err)
		}
		scatter.GlyphStyle.Shape = draw.CircleGlyph{}
		scatter.GlyphStyle.Radius = vg.Points(4)
		scatter.Color = theme.StatePoint()
		p.Add(scatter)

		if cfg.NumberStates {
			labels, _ := plotter.NewLabels(plotter.XYLabels{
				XYs:    []plotter.XY{{X: stateV, Y: state.Pressure}},
				Labels: []string{fmt.Sprintf("%d", i+1)},
			})
			labels.Offset.X = vg.Points(5)
			labels.Offset.Y = vg.Points(5)
			labels.TextStyle[0].Color = theme.StateNumber()

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
		width = DefaultPVConfig(cfg.Type).Width
	}
	height := cfg.Height
	if height == 0 {
		height = DefaultPVConfig(cfg.Type).Height
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
