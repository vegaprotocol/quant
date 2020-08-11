package pricedistribution

import (
	"math"
	"testing"

	"gonum.org/v1/gonum/stat/distuv"
)

func TestPriceRangeUniform(t *testing.T) {
	tolerance := 1e-12
	min := 100.0
	max := 200.0
	uniform := distuv.Uniform{Min: min, Max: max}
	alpha := 0.9
	expectedSmin := min + (max-min)*(1-alpha)/2
	expectedSmax := max - (max-min)*(1-alpha)/2

	Smin, Smax := PriceRange(uniform, alpha)

	alphaRecalculated := (Smax - Smin) / (max - min)
	if alphaRecalculated-alpha > tolerance {
		t.Logf("expected alpha=%g\n", alpha)
		t.Logf("actual alpha=%g\n", alphaRecalculated)
		t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(alphaRecalculated-alpha), tolerance)
	}
	if math.Abs(expectedSmin-Smin) > tolerance {
		t.Logf("expected S_min=%g\n", expectedSmin)
		t.Logf("actual S_min=%g\n", Smin)
		t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(expectedSmin-Smin), tolerance)
	}
	if math.Abs(expectedSmax-Smax) > tolerance {
		t.Logf("expected S_max=%g\n", expectedSmax)
		t.Logf("actual S_max=%g\n", Smax)
		t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(expectedSmax-Smax), tolerance)
	}
}

func TestPriceDistributionUniform(t *testing.T) {
	tolerance := 1e-12
	min := 100.0
	max := 200.0
	nBins := 1000
	increment := (max - min) / float64(nBins)
	bins := make([]float64, nBins+1)
	for i := range bins {
		bins[i] = min + float64(i)*increment
	}
	expectedProbabilityPerBin := 1 / float64(nBins)
	uniform := distuv.Uniform{Min: min, Max: max}

	probabilities := PriceDistribution(uniform, bins)

	totalProbability := 0.0
	for _, p := range probabilities {
		if math.Abs(expectedProbabilityPerBin-p) > tolerance {
			t.Logf("expected probability=%g\n", expectedProbabilityPerBin)
			t.Logf("actual probability=%g\n", p)
			t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(expectedProbabilityPerBin-p), tolerance)
		}
		totalProbability += p
	}

	if math.Abs(totalProbability-1) > tolerance {
		t.Log("expected total probability=1\n")
		t.Logf("actual total probability=%g\n", totalProbability)
		t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(totalProbability-1), tolerance)
	}
}
