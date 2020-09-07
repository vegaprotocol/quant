package riskmodelbs

import (
	"math"
	"testing"

	"code.vegaprotocol.io/quant/interfaces"
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
