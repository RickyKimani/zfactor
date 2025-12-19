package leekesler

import "testing"

func TestFindIndex(t *testing.T) {
	testCases := []struct {
		name string
		arr  []float64
		val  float64
		want int
	}{
		{
			name: "Standard case in the middle",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  35.0,
			want: 2,
		},
		{
			name: "Value is an exact lower bound of an interval",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  30.0,
			want: 1,
		},
		{
			name: "Value is an exact upper bound of an interval",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  40.0,
			want: 2,
		},
		{
			name: "Value in the very first interval",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  15.0,
			want: 0,
		},
		{
			name: "Value is the very first element",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  10.0,
			want: 0,
		},
		{
			name: "Value in the very last interval",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  45.0,
			want: 3,
		},
		{
			name: "Value is smaller than the entire range",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  5.0,
			want: 0,
		},
		{
			name: "Value is larger than the entire range",
			arr:  []float64{10, 20, 30, 40, 50},
			val:  55.0,
			want: 3,
		},
		{
			name: "Array with negative and zero values",
			arr:  []float64{-10, -5, 0, 5, 10},
			val:  -2.0,
			want: 1,
		},
		{
			name: "Empty slice",
			arr:  []float64{},
			val:  30.0,
			want: -1,
		},
		{
			name: "Slice with only one element",
			arr:  []float64{100},
			val:  100.0,
			want: -1,
		},
		{
			name: "Slice with two elements",
			arr:  []float64{10, 20},
			val:  15.0,
			want: 0,
		},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			got := findIndex(tc.arr, tc.val)
			if got != tc.want {
				t.Errorf("findIndex(%v, %f) = %d; want %d", tc.arr, tc.val, got, tc.want)
			}
		})
	}
}
