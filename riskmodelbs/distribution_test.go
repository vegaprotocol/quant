package riskmodelbs

import (
	"math"
	"testing"

	"code.vegaprotocol.io/quant/interfaces"
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

func TestLogNormalDistributionInternallyConsistent(t *testing.T) {
	const r float64 = 0.0
	const mu float64 = 0.2
	const sigma float64 = 1.5
	const S0 float64 = 123.0
	const tau = 1

	bsModelParameters := ModelParamsBS{Mu: 0, R: 0, Sigma: 1.5}

	distribution := bsModelParameters.GetProbabilityDistribution(S0, tau)

	ps := []float64{0, 0.001, 0.2, 0.5, 1, 10, 20, 30, 100, 1000, 10000}

	tolerance := bsModelParameters.GetProbabilityTolerance()
	for _, p := range ps {
		cdf := distribution.CDF(p)
		pRecalculated := distribution.Quantile(cdf)

		if math.Abs(p-pRecalculated) > tolerance {
			t.Errorf("Error=%g is more than tolerance (%g)", math.Abs(p-pRecalculated), tolerance)
		}
	}
}
