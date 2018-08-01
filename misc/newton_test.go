package misc

import (
	"math"
	"testing"
)

// TestNewtonViaSqrt we check whether newton converges to sqrt(2)
// by asking it to solve the nonlinear equation x^2 - 2 = 0;
func TestNewtonViaSqrt(t *testing.T) {
	const tolerance float64 = 1e-6

	f := func(x float64) float64 { return x*x - 2.0 }
	var x0 float64 = 1.5
	x, err := FindRootWithoutDerivative(f, x0, 100, 1e-7)
	if err != nil {
		t.Errorf("Newton solver failed\n")
	}

	error := math.Abs(math.Sqrt(2.0) - x)
	if math.IsNaN(error) || math.IsInf(error, 0) || error > tolerance {
		t.Errorf("Newton solver failed\n")
	}
}

func TestNewtonViaSqrtWithDerivative(t *testing.T) {
	const tolerance float64 = 1e-6

	f := func(x float64) float64 { return x*x - 2.0 }
	fPrime := func(x float64) float64 { return 2 * x }
	var x0 float64 = 1.5
	x, err := FindRoot(f, fPrime, x0, 100, 1e-7)
	if err != nil {
		t.Errorf("Newton solver failed\n")
	}

	error := math.Abs(math.Sqrt(2.0) - x)
	if math.IsNaN(error) || math.IsInf(error, 0) || error > tolerance {
		t.Errorf("Newton solver failed\n")
	}
}

func TestNewtonZeroDerivativeFail(t *testing.T) {
	f := func(x float64) float64 { return x*x - 2.0 }
	fPrime := func(x float64) float64 { return 2 * x }
	x, err := FindRoot(f, fPrime, 0.0, 100, 1e-7)
	if err == nil {
		t.Errorf("Newton solver shouldn't have worked in this case.\n")
		t.Logf("Returned x=%g", x)
	}
}

func TestNewtonNoConvergence(t *testing.T) {
	f := func(x float64) float64 { return x*x*x - 2*x + 2 }
	fPrime := func(x float64) float64 { return 3*x*x - 2 }
	x, err := FindRoot(f, fPrime, 1.0, 100, 1e-7)
	if err == nil {
		t.Errorf("Newton solver shouldn't have worked in this case.\n")
		t.Logf("Returned x=%g\n", x)
	}
	t.Logf("Error message is " + err.Error())
}
