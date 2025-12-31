package liquids

import (
	"testing"
)

func TestReducedDensity(t *testing.T) {
	tests := []struct {
		name    string
		tr      float64
		pr      float64
		wantErr bool
	}{
		// Valid Cases
		{
			name:    "Exact Match Tr=0.3 Pr=0.1",
			tr:      0.3,
			pr:      0.1,
			wantErr: false,
		},
		{
			name:    "Interpolation Tr=0.35 Pr=0.1",
			tr:      0.35,
			pr:      0.1,
			wantErr: false,
		},
		{
			name:    "Interpolation Tr=0.3 Pr=0.055",
			tr:      0.3,
			pr:      0.055,
			wantErr: false,
		},
		{
			name:    "Double Interpolation Tr=0.35 Pr=0.055",
			tr:      0.35,
			pr:      0.055,
			wantErr: false,
		},

		// Edge Cases - Temperature
		{
			name:    "Tr Too Low (0.2)",
			tr:      0.2,
			pr:      1.0,
			wantErr: true,
		},
		{
			name:    "Tr Too High (1.1)",
			tr:      1.1,
			pr:      1.0,
			wantErr: true,
		},

		// Edge Cases - Pressure
		{
			name:    "Pr Too High General (11.0)",
			tr:      0.5,
			pr:      11.0,
			wantErr: true,
		},
		{
			name:    "Pr High Fallback for Tr=0.97 (4.1)",
			tr:      0.97,
			pr:      4.1,
			wantErr: false,
		},
		{
			name:    "Pr Too Low (-0.1)",
			tr:      0.5,
			pr:      -0.1,
			wantErr: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ReducedDensity(tt.tr, tt.pr)
			if (err != nil) != tt.wantErr {
				t.Errorf("ReducedDensity() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !tt.wantErr && got <= 0 {
				t.Errorf("ReducedDensity() = %v, want > 0", got)
			}
		})
	}
}
