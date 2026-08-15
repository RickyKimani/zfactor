// Command phasor explores the "volume phasor": the complex-conjugate volume
// roots a cubic EOS produces when a state has only one real root.
//
// It maps beta = |Im V| over the CO2 (T, P) plane for the van der Waals EOS and
// overlays two independently-computed curves:
//
//   - the spinodal (dP/dV = 0), the mechanical stability limit, which is where
//     the cubic's discriminant changes sign and the complex pair is born; and
//   - the binodal (saturation pressure), the actual coexistence line.
//
// The thesis: beta -> 0 exactly on the spinodal (not the binodal), pinching to
// a point at the critical point. Between the binodal and the spinodal — the
// metastable region — beta is still zero, because three real roots survive
// there. So the phasor "turns on" at the limit of metastability.
//
// Run from the module root:  go run ./examples/phasor
package main

import (
	"fmt"
	"image/color"
	"log"
	"math"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/substance"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/palette"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
	"gonum.org/v1/plot/vg/draw"
)

// R in bar*cm^3/(mol*K), the unit convention used across zfactor examples.
const R = 10 * zfactor.RSI

// realTol is the |Im| below which a root counts as real.
const realTol = 1e-7

func main() {
	co2 := substance.CarbonDioxide
	Tc, Pc := co2.Critical.Tc, co2.Critical.Pc
	base := cubic.NewvdWCfg(0, 0, Tc, Pc, R)

	// phasor returns beta = |Im V| of the complex-conjugate root pair (0 when
	// all three roots are real) and how many of the three roots are real.
	phasor := func(T, P float64) (beta float64, nReal int) {
		cfg := *base
		cfg.T, cfg.P = T, P
		res, err := cubic.SolveForVolume(&cfg)
		if err != nil {
			return math.NaN(), 0
		}
		for _, v := range res.Volumes {
			im := math.Abs(imag(v))
			if im < realTol {
				nReal++
			}
			if im > beta {
				beta = im
			}
		}
		return beta, nReal
	}

	// --- the beta field over the (T, P) plane ---
	const nx, ny = 280, 280
	Tmin, Tmax := 220.0, 360.0
	Pmin, Pmax := 1.0, 140.0
	g := &grid{
		xs: linspace(Tmin, Tmax, nx),
		ys: linspace(Pmin, Pmax, ny),
		z:  make([][]float64, nx),
	}
	for i, T := range g.xs {
		g.z[i] = make([]float64, ny)
		for j, P := range g.ys {
			b, _ := phasor(T, P)
			g.z[i][j] = logFloor(b) // log scale so the near-transition gradient is visible
		}
	}

	// --- spinodal, exact vdW parametric trace: sweep V, get (T, P) ---
	// a, b match cubic.calculateA/B with alpha = 1 (vdW).
	a := (27.0 / 64.0) * R * R * Tc * Tc / Pc
	b := (1.0 / 8.0) * R * Tc / Pc
	var spin plotter.XYs
	for _, V := range linspace(b*1.0001, b*45, 1200) {
		T := 2 * a * (V - b) * (V - b) / (R * V * V * V) // dP/dV = 0  =>  T(V)
		P := R*T/(V-b) - a/(V*V)
		if T >= Tmin && T <= Tmax && P >= Pmin && P <= Pmax {
			spin = append(spin, plotter.XY{X: T, Y: P})
		}
	}

	// --- binodal, from the equal-fugacity saturation solver ---
	var bino plotter.XYs
	for _, T := range linspace(Tmin, Tc-0.5, 240) {
		cfg := *base
		Psat, err := cubic.SaturationPressure(&cfg, T)
		if err == nil && Psat >= Pmin && Psat <= Pmax {
			bino = append(bino, plotter.XY{X: T, Y: Psat})
		}
	}

	writePlot(g, spin, bino, Tc, Pc)
	printIsotherm(phasor, base, Tc, 280.0)
}

// writePlot renders the heat map with the spinodal / binodal overlays.
func writePlot(g *grid, spin, bino plotter.XYs, Tc, Pc float64) {
	p := plot.New()
	p.Title.Text = "CO2 volume phasor  |Im V|  (van der Waals, log scale)"
	p.X.Label.Text = "Temperature (K)"
	p.Y.Label.Text = "Pressure (bar)"

	p.Add(plotter.NewHeatMap(g, palette.Heat(24, 1)))

	spinLine, err := plotter.NewLine(spin)
	if err != nil {
		log.Fatal(err)
	}
	spinLine.Color = color.White
	spinLine.Width = vg.Points(2.2)

	binoLine, err := plotter.NewLine(bino)
	if err != nil {
		log.Fatal(err)
	}
	binoLine.Color = color.RGBA{R: 90, G: 200, B: 255, A: 255}
	binoLine.Width = vg.Points(2.2)
	binoLine.Dashes = []vg.Length{vg.Points(5), vg.Points(3)}

	cp, err := plotter.NewScatter(plotter.XYs{{X: Tc, Y: Pc}})
	if err != nil {
		log.Fatal(err)
	}
	cp.GlyphStyle.Color = color.White
	cp.GlyphStyle.Radius = vg.Points(4)
	cp.GlyphStyle.Shape = draw.CircleGlyph{}

	p.Add(spinLine, binoLine, cp)
	p.Legend.Add("spinodal dP/dV=0  (phasor -> 0)", spinLine)
	p.Legend.Add("binodal (saturation)", binoLine)
	p.Legend.Add("critical point", cp)
	p.Legend.Top = true

	const out = "examples/phasor/phasor_co2.png" // raster: a dense grid as SVG is ~7 MB of <rect>s
	if err := p.Save(22*vg.Centimeter, 16*vg.Centimeter, out); err != nil {
		log.Fatal(err)
	}
	fmt.Printf("wrote %s\n\n", out)
}

// printIsotherm scans P at fixed T, locating the spinodal pressures as the
// crossings where the real-root count changes, and marking the binodal.
func printIsotherm(phasor func(T, P float64) (float64, int), base *cubic.EOSCfg, Tc, Tiso float64) {
	cfg := *base
	Psat, _ := cubic.SaturationPressure(&cfg, Tiso)

	// Fine scan to find where nReal flips (the spinodal pressures).
	var spinodals []float64
	prev := 0
	for P := 1.0; P <= 130.0; P += 0.02 {
		_, n := phasor(Tiso, P)
		if prev != 0 && n != prev {
			spinodals = append(spinodals, P)
		}
		prev = n
	}

	fmt.Printf("CO2 van der Waals - phasor across the T = %.1f K isotherm (Tr = %.3f)\n", Tiso, Tiso/Tc)
	fmt.Printf("  binodal  (saturation)  P_sat = %.2f bar\n", Psat)
	for _, ps := range spinodals {
		fmt.Printf("  spinodal (nReal flips) P     = %.2f bar\n", ps)
	}
	fmt.Printf("\n  P (bar)   real roots   beta = |Im V| (cm^3/mol)\n")
	for P := 5.0; P <= 95.0; P += 10 {
		beta, n := phasor(Tiso, P)
		note := ""
		switch {
		case n == 3:
			note = "  three real roots -> phasor off"
		case P < Psat:
			note = "  outside lower spinodal (gas)"
		default:
			note = "  outside upper spinodal (liquid)"
		}
		fmt.Printf("  %6.1f      %d          %12.4f%s\n", P, n, beta, note)
	}
}

// grid implements plotter.GridXYZ: X is temperature, Y is pressure, Z is beta.
type grid struct {
	xs, ys []float64
	z      [][]float64
}

func (g *grid) Dims() (int, int)   { return len(g.xs), len(g.ys) }
func (g *grid) X(c int) float64    { return g.xs[c] }
func (g *grid) Y(r int) float64    { return g.ys[r] }
func (g *grid) Z(c, r int) float64 { return g.z[c][r] }

// logFloor maps beta to log10, flooring the three-real-root region (beta ~ 0)
// to a single background level so the color scale spans the live phasor field.
func logFloor(beta float64) float64 {
	const floor = -3.0
	if beta < 1e-3 {
		return floor
	}
	return math.Max(floor, math.Log10(beta))
}

// linspace returns n evenly spaced points on [lo, hi].
func linspace(lo, hi float64, n int) []float64 {
	out := make([]float64, n)
	if n == 1 {
		out[0] = lo
		return out
	}
	step := (hi - lo) / float64(n-1)
	for i := range out {
		out[i] = lo + step*float64(i)
	}
	return out
}
