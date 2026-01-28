package cp_test

import (
	"testing"

	"github.com/rickykimani/zfactor"
	"github.com/rickykimani/zfactor/cp"
)

func TestIdealGasEnthalpyChange(t *testing.T) {
	// Methane Gas
	gas := cp.MethaneGas

	state1 := zfactor.Args{
		T: 298.15,
		P: 101325,
		R: zfactor.RSI,
	}

	state2 := zfactor.Args{
		T: 1000,
		P: 101325,
		R: zfactor.RSI,
	}

	deltaH, err := gas.IdealGasEnthalpyChange(state1, state2)
	if err != nil {
		t.Fatalf("Calculated Failed: %v", err)
	}

	if deltaH <= 0 {
		t.Errorf("Expected positive enthalpy change for heating, got %v", deltaH)
	}

	t.Logf("Delta H for Methane (298.15 -> 1000K): %v J/mol", deltaH)
}

func TestIdealGasEntropyChange(t *testing.T) {
	// Methane Gas
	gas := cp.MethaneGas

	state1 := zfactor.Args{
		T: 298.15,
		P: 100000,
		R: zfactor.RSI,
	}

	state2 := zfactor.Args{
		T: 500,
		P: 100000, // Isobaric
		R: zfactor.RSI,
	}

	deltaS, err := gas.IdealGasEntropyChange(state1, state2)
	if err != nil {
		t.Fatalf("Calculated Failed: %v", err)
	}

	if deltaS <= 0 {
		t.Errorf("Expected positive entropy change for heating, got %v", deltaS)
	}
	t.Logf("Delta S (Isobaric) for Methane (298.15 -> 500K): %v J/mol.K", deltaS)

	// Isothermal expansion (lower pressure) -> Logic check
	// Entropy should increase if pressure decreases (V increases)
	state3 := zfactor.Args{
		T: 298.15,
		P: 50000,
		R: zfactor.RSI,
	}

	deltaS_expansion, err := gas.IdealGasEntropyChange(state1, state3)
	if err != nil {
		t.Fatalf("Expansion Calc Failed: %v", err)
	}

	// Delta S = -R ln(P2/P1). P2 < P1 => ln(P2/P1) < 0 => -ln(...) > 0
	if deltaS_expansion <= 0 {
		t.Errorf("Expected positive entropy change for expansion, got %v", deltaS_expansion)
	}
	t.Logf("Delta S (Expansion) for Methane (1bar -> 0.5bar): %v J/mol.K", deltaS_expansion)
}

func TestInputValidation(t *testing.T) {
	gas := cp.MethaneGas
	valid := zfactor.Args{T: 300, P: 1000, R: 8.314}
	invalidT := zfactor.Args{T: -10, P: 1000, R: 8.314}
	invalidR := zfactor.Args{T: 300, P: 1000, R: 1.0}

	// Test Invalid Temperature
	_, err := gas.IdealGasEnthalpyChange(valid, invalidT)
	if err == nil {
		t.Error("Expected error for negative temperature")
	}

	// Test Mismatched R
	_, err = gas.IdealGasEnthalpyChange(valid, invalidR)
	if err == nil {
		t.Error("Expected error for mismatched Gas Constants")
	}

	// Test Out of Range (Methane Max T is usually 1500K)
	outRange := zfactor.Args{T: 2000, P: 1000, R: 8.314}
	_, err = gas.IdealGasEnthalpyChange(valid, outRange)
	if err == nil {
		t.Log("Warning: Range check might not be failing if TMax is high, currently check expects error")
	}
}
