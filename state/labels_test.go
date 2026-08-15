package state_test

import (
	"os"
	"path/filepath"
	"regexp"
	"strconv"
	"testing"

	"github.com/rickykimani/zfactor/cubic"
	"github.com/rickykimani/zfactor/state"
	"github.com/rickykimani/zfactor/substance"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/font"
	"gonum.org/v1/plot/text"
	"gonum.org/v1/plot/vg"
)

// Isotherm labels are drawn just past the right end of their curve, which is
// where the x-axis ends. Unless the axis reserves room for them, the text
// starts at the edge of the plotting area and is cut off by the canvas
// boundary — the defect this guards against had "Tc=305.3 K" losing its
// second half.
//
// The check is made on the rendered SVG rather than on the axis limit,
// because the limit is only the means; whether the text lands on the canvas
// is the property worth holding. Widths come from the same font machinery
// the renderer used, so no glyph metrics are guessed here.
var (
	svgCanvasPattern = regexp.MustCompile(`width="(\d+)pt"`)
	svgTextPattern   = regexp.MustCompile(`(?s)<text x="([-0-9.]+)".*?font-size:([0-9.]+)px"[^>]*>([^<]*)</text>`)
	isothermPattern  = regexp.MustCompile(`^Tc?=[0-9.]+ K$`)
)

// isothermLabel is one label found in a rendered diagram.
type isothermLabel struct {
	text  string
	start float64 // leftmost point of the text, in points from the left edge
	width float64 // rendered width of the text, in points
}

// end returns the point the text finishes at.
func (l isothermLabel) end() float64 { return l.start + l.width }

// parseIsothermLabels returns the isotherm labels in an SVG, and the width of
// the canvas they have to fit inside.
func parseIsothermLabels(t *testing.T, path string) (labels []isothermLabel, canvas float64) {
	t.Helper()

	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("reading rendered diagram: %v", err)
	}

	source := string(raw)

	canvasMatch := svgCanvasPattern.FindStringSubmatch(source)
	if canvasMatch == nil {
		t.Fatal("no width attribute in the rendered SVG")
	}

	canvas, err = strconv.ParseFloat(canvasMatch[1], 64)
	if err != nil {
		t.Fatalf("parsing canvas width %q: %v", canvasMatch[1], err)
	}

	for _, match := range svgTextPattern.FindAllStringSubmatch(source, -1) {
		content := match[3]
		if !isothermPattern.MatchString(content) {
			continue
		}

		start, err := strconv.ParseFloat(match[1], 64)
		if err != nil {
			t.Fatalf("parsing label x %q: %v", match[1], err)
		}

		size, err := strconv.ParseFloat(match[2], 64)
		if err != nil {
			t.Fatalf("parsing font size %q: %v", match[2], err)
		}

		style := text.Style{
			Font:    font.From(plot.DefaultFont, vg.Length(size)),
			Handler: plot.DefaultTextHandler,
		}

		labels = append(labels, isothermLabel{
			text:  content,
			start: start,
			width: float64(style.Width(content)),
		})
	}

	return labels, canvas
}

// drawLabelledDiagram renders the two-state ethane diagram at one size, with
// isotherm labels on, and returns the path written.
func drawLabelledDiagram(t *testing.T, width, height state.Length) string {
	t.Helper()

	first, err := state.NewState(substance.Ethane, 299, 32)
	if err != nil {
		t.Fatalf("creating first state: %v", err)
	}

	second, err := state.NewState(substance.Ethane, 490, 70)
	if err != nil {
		t.Fatalf("creating second state: %v", err)
	}

	cfg := &state.PVConfig{
		Type:           &cubic.PR{},
		Title:          "PV Diagram for Ethane",
		LabelIsotherms: true,
		NumberStates:   true,
		Width:          width,
		Height:         height,
	}

	path := filepath.Join(t.TempDir(), "pv.svg")

	if err := state.DrawPV(cfg, path, first, second); err != nil {
		t.Fatalf("DrawPV: %v", err)
	}

	return path
}

func TestIsothermLabelsStayOnTheCanvas(t *testing.T) {
	// The default size and a range around it. The reserved room is a share
	// of the figure, so the small end is the demanding case: there the
	// label is a large fraction of the width.
	sizes := []struct {
		name          string
		width, height state.Length
	}{
		{"default 6x4", 6 * state.Inch, 4 * state.Inch},
		{"wide 10x6", 10 * state.Inch, 6 * state.Inch},
		{"wider 16x9", 16 * state.Inch, 9 * state.Inch},
	}

	for _, size := range sizes {
		t.Run(size.name, func(t *testing.T) {
			labels, canvas := parseIsothermLabels(t, drawLabelledDiagram(t, size.width, size.height))

			// Both isotherms and the critical isotherm are labelled, so
			// finding fewer than three means the labels stopped being
			// drawn or the parsing stopped matching. Either way the
			// assertions below would pass vacuously.
			const wantLabels = 3

			if len(labels) != wantLabels {
				t.Fatalf("found %d isotherm labels, want %d: the diagram or the parsing changed",
					len(labels), wantLabels)
			}

			// A real gap is required rather than merely "not past the
			// edge". Without the reservation these labels land within
			// hundredths of a point of the boundary, which satisfies
			// end <= canvas yet renders as clipped text, so that
			// weaker assertion would not notice the defect at all.
			//
			// The reservation currently leaves at least 5 pt at this
			// size, so the threshold has room to spare while staying
			// well clear of the flush case.
			const minClearance = 2.0

			for _, label := range labels {
				if gap := canvas - label.end(); gap < minClearance {
					t.Errorf("label %q has %.2f pt of clearance, want at least %.2f: starts at %.2f, is %.2f wide, ends at %.2f on a %.2f canvas",
						label.text, gap, minClearance, label.start, label.width, label.end(), canvas)
				}
			}
		})
	}
}

// The room reserved is a share of the figure width, so a larger figure should
// leave the labels further clear of the edge. This pins the scaling rather
// than any single measurement, which is what keeps the fix working at sizes
// the test above does not name.
func TestIsothermLabelMarginGrowsWithFigureWidth(t *testing.T) {
	margin := func(width, height state.Length) float64 {
		labels, canvas := parseIsothermLabels(t, drawLabelledDiagram(t, width, height))
		if len(labels) == 0 {
			t.Fatal("no isotherm labels found")
		}

		// The tightest label is the one that decides whether anything is
		// clipped.
		worst := canvas - labels[0].end()
		for _, label := range labels[1:] {
			if m := canvas - label.end(); m < worst {
				worst = m
			}
		}

		return worst
	}

	small := margin(6*state.Inch, 4*state.Inch)
	large := margin(16*state.Inch, 9*state.Inch)

	if small <= 0 {
		t.Errorf("no margin at the default size: %.2f pt", small)
	}

	if large <= small {
		t.Errorf("margin did not grow with the figure: %.2f pt at 6in, %.2f pt at 16in", small, large)
	}
}
