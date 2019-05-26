package riskmodelbsjmp

import (
	"fmt"
	"math"
	"testing"
	"time"

	"gonum.org/v1/gonum/stat/distuv"
)

const testTolerance float64 = 1.0e-8

// TestBSFwdRiskFactorsVal is a most basic test where we compare to a value we
// obtained from IPython.
func TestBSFwdRiskFactorsVal(t *testing.T) {
	const r float64 = 0.0
	const mu float64 = 0.0
	const sigma float64 = 1.5

	const tau = 1.0 / 365 / 24 / 60
	const lambda float64 = 0.01

	var paramsBs RiskModelParamsBSJmp
	paramsBs.mu = mu
	paramsBs.r = r
	paramsBs.sigma = sigma
	riskFactors := RiskFactorsForward(lambda, tau, paramsBs)
	const rfShort float64 = 5.5276558362e-03 // from IPython notebook
	const rfLong float64 = 5.5011137953e-03

	error := math.Abs(riskFactors.Short-rfShort) + math.Abs(riskFactors.Long-rfLong)
	if math.IsNaN(error) || math.IsInf(error, 0) || error > testTolerance {
		t.Errorf("Error %g is bigger than tolerance %g\n", error, testTolerance)
	}
}

// TestBSFwdRiskFactorsOrder makes sure that the risk factor for long
// positions doesn't exceed the risk factor for short positions.
func TestBSFwdRiskFactorsOrder(t *testing.T) {
	const lambda float64 = 0.01
	const numRuns int = 1000

	for runIdx := 0; runIdx < numRuns; runIdx++ {
		var paramsBs RiskModelParamsBSJmp
		paramsBs.mu = distuv.UnitUniform.Rand()
		paramsBs.r = distuv.UnitUniform.Rand()
		paramsBs.sigma = distuv.UnitUniform.Rand()
		tau := distuv.UnitUniform.Rand()

		riskFactors := RiskFactorsForward(lambda, tau, paramsBs)

		if math.IsNaN(riskFactors.Short) || math.IsNaN(riskFactors.Long) ||
			math.IsInf(riskFactors.Short, 0) || math.IsInf(riskFactors.Long, 0) {
			t.Errorf("risk factor is NaN or Inf")
		}

		if (riskFactors.Short + testTolerance) < riskFactors.Long {
			t.Logf("r=%g, mu=%g, sigma=%g, tau=%g\n", paramsBs.r, paramsBs.mu, paramsBs.sigma, tau)
			t.Logf("rf short=%g, rf long=%g\n", riskFactors.Short, riskFactors.Long)
			t.Errorf("risk factor for long must be less than that for short")
		}
	}
}

func TestTimeTakenForForwardsRiskFactor(t *testing.T) {
	var paramsBs RiskModelParamsBSJmp
	paramsBs.mu = 0.1
	paramsBs.r = 0.01
	paramsBs.sigma = 1.2

	const numRuns int = 1000000
	start := time.Now()

	for i := 0; i < numRuns; i++ {
		RiskFactorsForward(0.01, 1.0/365.25/24, paramsBs)
	}

	elapsed := time.Since(start)

	fmt.Printf("Num of times we can calculate forward risk factors in BS model: %.1f per second.\n",
		float64(numRuns)/elapsed.Seconds())
}
