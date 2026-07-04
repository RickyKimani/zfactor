package state

import "errors"

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
