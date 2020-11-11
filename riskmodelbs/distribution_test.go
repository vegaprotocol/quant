package riskmodelbs

import (
	"math"
	"sort"
	"testing"

	"code.vegaprotocol.io/quant/interfaces"
	"code.vegaprotocol.io/quant/pricedistribution"
	"gonum.org/v1/gonum/stat"
	"gonum.org/v1/gonum/stat/distuv"
)

func TestBSIsAnalyticalModel(t *testing.T) {
	const r float64 = 0.0
	const mu float64 = 0.0
	const sigma float64 = 1.5

	var analyticBs interfaces.AnalyticalModel = ModelParamsBS{Mu: 0, R: 0, Sigma: 1.5}

	if analyticBs == nil {
		t.Error("Expeced ModelParamsBS to implement AnalyticalModel interface")
	}
}

func TestLogNormalMeanAndVariace(t *testing.T) {
	tolerance := 1e-6
	r := 0.0
	S := 10.0
	mu := 0.1
	tau := 0.25
	sigma := 1.0
	expectedMean := S * math.Exp(mu*tau)
	expectedVariance := S * S * math.Exp(2*mu*tau) * (math.Exp(sigma*sigma*tau) - 1)

	dist := ModelParamsBS{Mu: mu, R: r, Sigma: sigma}.GetProbabilityDistribution(S, tau)

	actualMean := dist.Mean()
	actualVariance := dist.Variance()

	if math.Abs(expectedMean-actualMean) > tolerance {
		t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(expectedMean-actualMean), tolerance)
	}

	if math.Abs(expectedVariance-actualVariance) > tolerance {
		t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(expectedVariance-actualVariance), tolerance)
	}
}

func TestLogNormalDistributionInternallyConsistent(t *testing.T) {
	const r float64 = 0.0
	const mu float64 = 0.2
	const sigma float64 = 1.5
	const S0 float64 = 123.0
	const tau = 1

	bsModelParameters := ModelParamsBS{Mu: mu, R: r, Sigma: sigma}

	distribution := bsModelParameters.GetProbabilityDistribution(S0, tau)

	xs := []float64{0, 0.001, 0.2, 0.5, 1, 10, 20, 30, 100, 1000, 10000}

	tolerance := bsModelParameters.GetProbabilityTolerance()
	for _, x := range xs {
		cdf := distribution.CDF(x)
		xRecalculated := distribution.Quantile(cdf)

		if math.Abs(x-xRecalculated) > tolerance {
			t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(x-xRecalculated), tolerance)
		}
	}
}

func TestLogNormalDistributionAgainstNormal(t *testing.T) {
	relativeTolerance := 1e-6
	const r float64 = 0.0
	const mu float64 = 0.0
	const sigma float64 = 2
	const S0 float64 = 110.0
	const tau = 1.0 / 60.0 / 24.0 / 365.25
	xs := []float64{0, 0.001, 0.1, 1, 10, 15, 23, 51, 100, 109, 110, 111, 163.5, 200, 500, 1000}
	bsModelParameters := ModelParamsBS{Mu: mu, R: r, Sigma: sigma}
	stdNormal := distuv.UnitNormal

	distribution := bsModelParameters.GetProbabilityDistribution(S0, tau)

	for _, x := range xs {
		y := (math.Log(x/S0) - (mu-0.5*sigma*sigma)*tau) / (sigma * math.Sqrt(tau))
		logNormalCdfX := distribution.CDF(x)
		stdNormalCdfY := stdNormal.CDF(y)

		if math.Abs(logNormalCdfX/stdNormalCdfY-1) > relativeTolerance {
			t.Errorf("Error=%g is more than tolerance (%g)(x=%g)", math.Abs(logNormalCdfX/stdNormalCdfY-1), relativeTolerance, x)
		}
	}
}

func Test_PriceRangeConsistentWithPriceDistribution(t *testing.T) {
	tolerance := 1e-6
	const r float64 = 0.016
	const mu float64 = 0.0
	const sigma float64 = 2
	const S0 float64 = 100000

	bsModelParameters := ModelParamsBS{Mu: mu, R: r, Sigma: sigma}

	const horizon = 1.0 / 24.0 / 365.25
	p := 0.9999

	distribution := bsModelParameters.GetProbabilityDistribution(S0, horizon)
	min, max := pricedistribution.PriceRange(distribution, p)

	probs := pricedistribution.PriceDistribution(distribution, []float64{min, max})

	if len(probs) != 1 {
		t.Error("Expected a single entry in the result array")
	}

	if math.Abs(p-probs[0]) > tolerance {
		t.Errorf("Error=%g is more than tolerance (%g)(p=%g)", math.Abs(p-probs[0]), tolerance, p)
	}
}

func GenerateAntitheticSamples(numIndepMCSamples int, bsModelParameters ModelParamsBS, S0 float64, tau float64) (probabilities []float64) {
	mu := bsModelParameters.Mu
	sigma := bsModelParameters.Sigma

	numMCSamples := 2 * numIndepMCSamples
	SatTau := make([]float64, numMCSamples)
	for i := 0; i < numIndepMCSamples; i++ {
		z := distuv.UnitNormal.Rand()
		SatTau[i] = S0 * math.Exp((mu-0.5*sigma*sigma)*tau+sigma*math.Sqrt(tau)*z)
		SatTau[numIndepMCSamples+i] = S0 * math.Exp((mu-0.5*sigma*sigma)*tau+sigma*math.Sqrt(tau)*(-z)) // antithetic sample
	}
	return SatTau
}

func TestLogNormalDistributionAgainstMonteCarlo(t *testing.T) {
	tolerance := 1e-3
	const numIndepMCSamples int = 2000000

	const r float64 = 0.0
	const mu float64 = 0.0
	const sigma float64 = 2
	const S0 float64 = 110.0
	const tau = 1.0 / 365.25
	xs := []float64{0, 0.001, 0.1, 1, 10, 15, 23, 51, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 163.5, 200, 500, 1000}

	bsModelParameters := ModelParamsBS{Mu: mu, R: r, Sigma: sigma}
	distribution := bsModelParameters.GetProbabilityDistribution(S0, tau)
	analyticProbabilities := pricedistribution.PriceDistribution(distribution, xs)
	// generate samples from lognormal distribution and sort them for use by CDF
	SatTau := GenerateAntitheticSamples(numIndepMCSamples, bsModelParameters, S0, tau)
	sort.Float64s(SatTau)

	for i := 0; i < len(xs)-1; i++ {
		cdfOfXLb := stat.CDF(xs[i], 1, SatTau, nil)
		cdfOfXUb := stat.CDF(xs[i+1], 1, SatTau, nil)
		probOfBin := cdfOfXUb - cdfOfXLb
		analyticProbOfBin := analyticProbabilities[i]
		diff := math.Abs(analyticProbOfBin - probOfBin)
		if diff > tolerance {
			t.Errorf("xs=%g, empirical prob of bin=%g, analytic prob of bin=%g, diff=%g \n", xs[i], probOfBin, analyticProbOfBin, diff)
		}
	}
}

func TestProbOfTradingAgainstMonteCarlo(t *testing.T) {
	relativeTolerance := 1e-3
	const numIndepMCSamples int = 4000000

	const r float64 = 0.0
	const mu float64 = 0.0
	const sigma float64 = 2
	const S0 float64 = 110.0
	const tau = 1.0 / 365.25
	xs := []float64{106, 107, 108, 109, 109.5, 110, 110.5, 111, 112, 113, 114, 115, 163.5, 200, 500, 1000}

	bsModelParameters := ModelParamsBS{Mu: mu, R: r, Sigma: sigma}
	distribution := bsModelParameters.GetProbabilityDistribution(S0, tau)

	// generate samples from lognormal distribution and sort them for use by CDF
	SatTau := GenerateAntitheticSamples(numIndepMCSamples, bsModelParameters, S0, tau)
	sort.Float64s(SatTau)

	for i := 0; i < len(xs); i++ {
		isBuySide := xs[i] <= S0
		analProb := pricedistribution.ProbabilityOfTrading(distribution, xs[i], isBuySide, false, 0, math.Inf(1))
		empProb := 0.0
		if isBuySide {
			empProb = stat.CDF(xs[i], 1, SatTau, nil)
		} else {
			empProb = 1.0 - stat.CDF(xs[i], 1, SatTau, nil)
		}
		diff := math.Abs(empProb - analProb)
		if diff > relativeTolerance {
			t.Errorf("xs=%g, empirical prob=%g, analytic prob of bin=%g, diff=%g \n", xs[i], empProb, analProb, diff)
		}
	}
}

func TestPriceRangeAgainstMonteCarlo(t *testing.T) {
	relativeTolerance := 1e-3
	const numIndepMCSamples int = 4000000

	const r float64 = 0.0
	const mu float64 = 0.0
	const sigma float64 = 2
	const S0 float64 = 12345.0

	taus := []float64{1.0 / 365.25, 1.0 / 365.25 / 24, 1.0 / 365.25 / 24 / 60, 1.0 / 365.25 / 24 / 60 / 60}
	alphas := []float64{0.5, 0.75, 0.875, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999}

	bsModelParameters := ModelParamsBS{Mu: mu, R: r, Sigma: sigma}

	for _, tau := range taus {
		distribution := bsModelParameters.GetProbabilityDistribution(S0, tau)
		// generate samples at tau
		SatTau := GenerateAntitheticSamples(numIndepMCSamples, bsModelParameters, S0, tau)
		sort.Float64s(SatTau)

		for _, alpha := range alphas {

			min, max := pricedistribution.PriceRange(distribution, alpha)
			empAlpha := getEmpricalAlphaFromRange(SatTau, min, max)
			relDiff := math.Abs(alpha/empAlpha - 1)
			if relDiff > relativeTolerance {
				t.Errorf("empirical alpha for range [%v, %v] is %v, while %v was prescribed, diff=%g \n", min, max, empAlpha, alpha, relDiff)
			}
		}
	}
}

func getEmpricalAlphaFromRange(allPrices []float64, min, max float64) float64 {
	pricesInRange := getPricesInRange(allPrices, min, max)

	return float64(len(pricesInRange)) / float64(len(allPrices))
}

func getPricesInRange(prices []float64, min, max float64) []float64 {
	iMin := 0
	iMax := len(prices) - 1
	for i := 0; i < len(prices); i++ {
		if prices[i] >= min {
			iMin = i
			break
		}
	}

	for i := iMax; i > iMin; i-- {
		if prices[i] <= max {
			iMax = i
			break
		}
	}

	return prices[iMin:iMax]
}
