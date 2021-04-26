package riskmeasures

import (
	"sort"

	"gonum.org/v1/gonum/stat"
)

//EmpiricalVaR Calculates empirical value at risk
func EmpiricalVaR(x []float64, alpha float64, isSorted bool) float64 {
	if !isSorted {
		sort.Float64s(x)
	}
	return -stat.Quantile(alpha, 1, x, nil)
}

// EmpiricalEs calculates the empirical expected shortfall for samples x
// It assumes x has already been sorted
func EmpiricalEs(x []float64, lambda float64, isSorted bool) float64 {
	N := len(x)
	empVar := EmpiricalVaR(x, lambda, isSorted)

	var countSamplesLessThanMinusVar int   // = 0 according to editor hints
	var sumSamplesLessThanMinusVar float64 // = 0.0 according to editor hints
	for i := 0; i < N; i++ {
		if x[i] <= -empVar {
			countSamplesLessThanMinusVar++
			sumSamplesLessThanMinusVar += x[i]
		}

	}
	pXleqMinusVar := float64(countSamplesLessThanMinusVar) / float64(N)
	es := (-1.0 / lambda) * (sumSamplesLessThanMinusVar/float64(N) + empVar*(pXleqMinusVar-lambda))
	//sigmaSq = stat.Variance(z, nil)
	return es
}
