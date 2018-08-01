package misc

import (
	"math"
	"testing"

	"gonum.org/v1/gonum/integrate"
	"gonum.org/v1/gonum/stat/distuv"
)

// Gonum documentation https://godoc.org/gonum.org/v1/gonum

func precomputeDensityValues(minX, maxX float64, numPoints int) ([]float64, []float64) {
	h := (maxX - minX) / float64(numPoints)
	x := make([]float64, numPoints)
	f := make([]float64, numPoints)
	for i := 0; i < numPoints; i++ {
		x[i] = minX + float64(i)*h
		f[i] = GaussDensity(x[i])
	}
	return x, f
}

func TestDensityViaIntegration(t *testing.T) {
	const tolerance float64 = 1e-10
	// these define how finely we do the numerical intergration
	const numIntegrationPoints int = 10000
	const minValXIntegration float64 = -7.0
	const maxValXIntegration float64 = 7

	x, f := precomputeDensityValues(minValXIntegration, maxValXIntegration, numIntegrationPoints)
	p := integrate.Trapezoidal(x, f)
	error := math.Abs(p - 1.0)
	if math.IsNaN(error) || math.IsInf(error, 0) || error > tolerance {
		t.Errorf("Gaussian density not intergating to 1 over the real line.\n")
	}
}

func TestCDFApproxViaGonum(t *testing.T) {
	const tolerance float64 = 1e-4
	// these define where we compare the CDF
	const numPoints int = 10000
	const minX float64 = -7.0
	const maxX float64 = 7
	h := (maxX - minX) / float64(numPoints)
	for i := 0; i < numPoints; i++ {
		x := minX + float64(i)*h
		cfdApprox := ApproxGaussCdf(x)
		cdfGonum := distuv.UnitNormal.CDF(x)
		error := math.Abs(cfdApprox - cdfGonum)
		if math.IsNaN(error) || math.IsInf(error, 0) || error > tolerance {
			t.Errorf("Gaussian CDF approx at x=%g is too far from Gonum value.\n", x)
		}
	}
}
