package abbott

import (
	"math"
	"testing"

	"github.com/rickykimani/zfactor"
)

func TestB0(t *testing.T) {
	tests := []struct {
		name    string
		tr      float64
		want    float64
		wantErr error
	}{
		{"Valid Tr=1", 1.0, 0.083 - 0.422, nil},
		{"Valid Tr=2", 2.0, 0.083 - 0.422/math.Pow(2, 1.6), nil},
		{"Invalid Tr=0", 0.0, 0, zfactor.ErrInvalidTr},
		{"Invalid Tr=-1", -1.0, 0, zfactor.ErrInvalidTr},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := B0(tt.tr)
			if err != tt.wantErr {
				t.Errorf("B0() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err == nil && math.Abs(got-tt.want) > 1e-9 {
				t.Errorf("B0() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestB1(t *testing.T) {
	tests := []struct {
		name    string
		tr      float64
		want    float64
		wantErr error
	}{
		{"Valid Tr=1", 1.0, 0.139 - 0.172, nil},
		{"Valid Tr=2", 2.0, 0.139 - 0.172/math.Pow(2, 4.2), nil},
		{"Invalid Tr=0", 0.0, 0, zfactor.ErrInvalidTr},
		{"Invalid Tr=-1", -1.0, 0, zfactor.ErrInvalidTr},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := B1(tt.tr)
			if err != tt.wantErr {
				t.Errorf("B1() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err == nil && math.Abs(got-tt.want) > 1e-9 {
				t.Errorf("B1() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestDB0(t *testing.T) {
	tests := []struct {
		name    string
		tr      float64
		want    float64
		wantErr error
	}{
		{"Valid Tr=1", 1.0, 0.675, nil},
		{"Valid Tr=2", 2.0, 0.675 / math.Pow(2, 2.6), nil},
		{"Invalid Tr=0", 0.0, 0, zfactor.ErrInvalidTr},
		{"Invalid Tr=-1", -1.0, 0, zfactor.ErrInvalidTr},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := DB0(tt.tr)
			if err != tt.wantErr {
				t.Errorf("DB0() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err == nil && math.Abs(got-tt.want) > 1e-9 {
				t.Errorf("DB0() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestDB1(t *testing.T) {
	tests := []struct {
		name    string
		tr      float64
		want    float64
		wantErr error
	}{
		{"Valid Tr=1", 1.0, 0.722, nil},
		{"Valid Tr=2", 2.0, 0.722 / math.Pow(2, 5.2), nil},
		{"Invalid Tr=0", 0.0, 0, zfactor.ErrInvalidTr},
		{"Invalid Tr=-1", -1.0, 0, zfactor.ErrInvalidTr},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := DB1(tt.tr)
			if err != tt.wantErr {
				t.Errorf("DB1() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if err == nil && math.Abs(got-tt.want) > 1e-9 {
				t.Errorf("DB1() = %v, want %v", got, tt.want)
			}
		})
	}
}
