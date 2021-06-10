package pricedistribution

import (
	"math"
	"testing"

	"code.vegaprotocol.io/quant/riskmodelbs"
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

	assert(t, "S_min", expectedSmin, Smin, tolerance)
	assert(t, "S_max", expectedSmax, Smax, tolerance)
	assert(t, "alpha", alpha, alphaRecalculated, tolerance)
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
		assert(t, "probability", expectedProbabilityPerBin, p, tolerance)
		totalProbability += p
	}

	assert(t, "total probability", 1, totalProbability, tolerance)
}

func TestProbabilityOfTradingUniformWithMinMax(t *testing.T) {
	tolerance := 1e-12
	lb := 100.0
	ub := 200.0
	uniform := distuv.Uniform{Min: lb, Max: ub}
	Smin := lb + 10
	Smax := ub - 10
	maxProb := 1.0

	var testCases = []struct {
		price    float64
		isBid    bool
		expected float64
	}{
		{lb, true, 0},
		{lb, false, 0},
		{Smin, true, 0},
		{Smin, false, maxProb},
		{Smin + tolerance, false, maxProb},
		{(Smax + Smin) / 2, false, maxProb / 2},
		{(Smax + Smin) / 2, true, maxProb / 2},
		{Smax, true, maxProb},
		{Smax - tolerance, true, maxProb},
		{Smax, false, 0},
		{ub, true, 0},
		{ub, false, 0},
	}

	for _, c := range testCases {
		actual := ProbabilityOfTrading(uniform, c.price, c.isBid, true, Smin, Smax)
		assert(t, "probability of trading", c.expected, actual, tolerance)
	}
}

func TestProbabilityOfTradingUniformNoMinMax(t *testing.T) {
	tolerance := 1e-12
	lb := 100.0
	ub := 200.0
	uniform := distuv.Uniform{Min: lb, Max: ub}

	var testCases = []struct {
		price    float64
		isBid    bool
		expected float64
	}{
		{lb, true, 0},
		{lb, false, 1},
		{(lb + ub) / 2, false, 0.5},
		{(lb + ub) / 2, true, 0.5},
		{ub, true, 1},
		{ub, false, 0},
	}

	for _, c := range testCases {
		actual := ProbabilityOfTrading(uniform, c.price, c.isBid, false, 0, 0)
		assert(t, "probability of trading", c.expected, actual, tolerance)
	}
}

func TestProbabilityOfTradingNormalisation(t *testing.T) {
	tolernace := 1e-2
	s0 := 100.0
	tau := 1 / 365.25 / 24
	min := 95.0
	max := 105.0
	sBid := 99.999
	sAsk := 100.001
	expectedProb := 0.5

	bsModel := riskmodelbs.ModelParamsBS{Mu: 0, R: 0, Sigma: 1.2}
	pdf := bsModel.GetProbabilityDistribution(s0, tau)

	prob1Bid := ProbabilityOfTrading(pdf, sBid, true, true, min, max)
	prob2Bid := ProbabilityOfTrading(pdf, sBid, true, false, math.NaN(), math.NaN())
	prob1Ask := ProbabilityOfTrading(pdf, sAsk, false, true, min, max)
	prob2Ask := ProbabilityOfTrading(pdf, sAsk, false, false, math.NaN(), math.NaN())

	assert(t, "probability of trading", expectedProb, prob1Bid, tolernace)
	assert(t, "probability of trading", expectedProb, prob2Bid, tolernace)
	assert(t, "probability of trading", prob1Bid, prob2Bid, tolernace)
	assert(t, "probability of trading", expectedProb, prob1Ask, tolernace)
	assert(t, "probability of trading", expectedProb, prob2Ask, tolernace)
	assert(t, "probability of trading", prob1Ask, prob2Ask, tolernace)

}

func assert(t *testing.T, label string, expected, actual, tolerance float64) {
	if math.Abs(expected-actual) > tolerance {
		t.Logf("expected %s=%g\n", label, expected)
		t.Logf("actual %s=%g\n", label, actual)
		t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(expected-actual), tolerance)
	}
}
