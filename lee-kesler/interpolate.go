package leekesler

import (
	"errors"
	"sort"
)

// At returns the interpolated value at the given reduced pressure (pr)
// and reduced temperature (tr). Returns an error if pr or tr are out of range.
//
// Usage:
//
//	v, err := leekesler.Z0Table.At(0.66, 1.2)
func (t table) At(pr, tr float64) (float64, error) {
	return interpolate(pr, tr, t)
}

// interpolate performs bilinear interpolation on the provided table.
// Returns an error if pr or tr are out of range.
func interpolate(pr, tr float64, table table) (float64, error) {
	//bounds
	if pr < table.Pr[0] || pr > table.Pr[len(table.Pr)-1] {
		return 0, errors.New("reduced pressure out of range")
	}
	if tr < table.Tr[0] || tr > table.Tr[len(table.Tr)-1] {
		return 0, errors.New("reduced temperature out of range")
	}

	i := findIndex(table.Pr, pr)
	j := findIndex(table.Tr, tr)

	x1, x2 := table.Pr[i], table.Pr[i+1]
	y1, y2 := table.Tr[j], table.Tr[j+1]

	// Values is organized as Values[TrIndex][PrIndex]
	M11 := table.Values[j][i]
	M12 := table.Values[j][i+1]
	M21 := table.Values[j+1][i]
	M22 := table.Values[j+1][i+1]

	M1 := ((x2-pr)/(x2-x1))*M11 + ((pr-x1)/(x2-x1))*M12
	M2 := ((x2-pr)/(x2-x1))*M21 + ((pr-x1)/(x2-x1))*M22

	M := ((y2-tr)/(y2-y1))*M1 + ((tr-y1)/(y2-y1))*M2

	return M, nil

}

func findIndex(arr []float64, val float64) int {
	if len(arr) < 2 {
		return -1
	}
	i := sort.SearchFloat64s(arr, val)

	if i == len(arr) {
		return len(arr) - 2
	}
	if i == 0 {
		return 0
	}

	return i - 1
}
