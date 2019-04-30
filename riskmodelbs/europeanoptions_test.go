package riskmodelbs

import (
	"fmt"
	"math"
	"testing"
	"time"

	"code.vegaprotocol.io/quant/bsformula"
	"code.vegaprotocol.io/quant/riskmeasures"

	"gonum.org/v1/gonum/stat/distuv"
)

// the inputs are S, K, mu, sigma, T, lambda, tau
var testValues = []struct {
	S      float64
	K      float64
	mu     float64
	sigma  float64
	T      float64
	lambda float64
	tau    float64
}{
	{1.0, 1.0, 0.05, 0.5, 0.25, 0.01, 1.0 / 365.25 / 24 / 60}, // test at the money
	{2.0, 1.0, -0.1, 0.5, 0.25, 0.01, 1.0 / 365.25 / 24 / 60},
	{1.0, 1.5, 0.07, 0.5, 0.25, 0.01, 1.0 / 365.25 / 24 / 60}, // test high strike
	{1.0, 0.5, 0.12, 0.5, 0.25, 0.01, 1.0 / 365.25 / 24 / 60}, // test low strike
	{1.0, 0.5, -0.2, 0.5, 0.25, 0.01, 1.0 / 365.25 / 24 / 60}, // test zero interest rate
	{1.0, 0.8, 0.00, 0.1, 2.00, 0.01, 1.0 / 365.25 / 24 / 60}, // test long maturity
	{75., 80., 0.03, 0.1, 2.00, 0.01, 1.0 / 365.25 / 24 / 60}, // test non-unit strike / mat
}

func TestCallPutRiskFactorsUsingMonteCarlo(t *testing.T) {
	const testToleranceForMC float64 = 1.0e-3
	const numIndepMCSamples int = 20000

	// generate normal samples for use later, using antithetic sampling
	numMCSamples := 2 * numIndepMCSamples
	Z := make([]float64, numMCSamples)
	for i := 0; i < numIndepMCSamples; i++ {
		z := distuv.UnitNormal.Rand()
		Z[i] = z
		Z[numIndepMCSamples+i] = -z // antithetic sample
	}

	for _, vals := range testValues {
		S := vals.S
		K := vals.K
		r := 0.0
		mu := vals.mu
		sigma := vals.sigma
		T := vals.T
		lambda := vals.lambda
		tau := vals.tau

		simLongCalls := make([]float64, numMCSamples)
		simShortCalls := make([]float64, numMCSamples)
		simLongPuts := make([]float64, numMCSamples)
		simShortPuts := make([]float64, numMCSamples)
		currentCall := bsformula.BSCallPrice(S, K, r, sigma, T)
		currentPut := bsformula.BSPutPrice(S, K, r, sigma, T)
		for i := 0; i < numMCSamples; i++ {
			SatTau := S * math.Exp((mu-0.5*sigma*sigma)*tau+sigma*math.Sqrt(tau)*Z[i])
			simLongCalls[i] = bsformula.BSCallPrice(SatTau, K, r, sigma, T-tau) - currentCall
			simShortCalls[i] = -simLongCalls[i]
			simLongPuts[i] = bsformula.BSPutPrice(SatTau, K, r, sigma, T-tau) - currentPut
			simShortPuts[i] = -simLongPuts[i]
		}
		empMarginLongCall := riskmeasures.EmpiricalEs(simLongCalls, lambda, false)
		empMarginShortCall := riskmeasures.EmpiricalEs(simShortCalls, lambda, false)
		empMarginLongPut := riskmeasures.EmpiricalEs(simLongPuts, lambda, false)
		empMarginShortPut := riskmeasures.EmpiricalEs(simShortPuts, lambda, false)

		riskFactorsCall := RiskFactorsCall(lambda, tau, S, K, T, ModelParamsBS{mu, r, sigma})

		marginShortCall := S * riskFactorsCall.Short
		marginLongCall := S * riskFactorsCall.Long

		errorCall := math.Abs(marginShortCall-empMarginShortCall) + math.Abs(marginLongCall-empMarginLongCall)
		errorCall /= currentCall
		if math.IsNaN(errorCall) || math.IsInf(errorCall, 0) || errorCall > testToleranceForMC {
			t.Logf("mu=%g, sigma=%g, tau=%g\n", mu, sigma, tau)
			t.Logf("margin short call=%g, emp margin short call=%g\n", marginShortCall, empMarginShortCall)
			t.Logf("margin long call=%g, emp margin long call=%g\n", marginLongCall, empMarginLongCall)
			t.Errorf("Error=%g is more than tolerance", errorCall)
		}

		riskFactorsPut := RiskFactorsPut(lambda, tau, S, K, T, ModelParamsBS{mu, r, sigma})

		marginShortPut := S * riskFactorsPut.Short
		marginLongPut := S * riskFactorsPut.Long

		errorPut := math.Abs(marginShortPut-empMarginShortPut) + math.Abs(marginLongPut-empMarginLongPut)
		errorPut /= currentPut
		if math.IsNaN(errorPut) || math.IsInf(errorPut, 0) || errorPut > testToleranceForMC {
			t.Logf("mu=%g, sigma=%g, tau=%g\n", mu, sigma, tau)
			t.Logf("margin short put=%g, emp margin long put=%g\n", marginShortPut, empMarginShortPut)
			t.Logf("margin long put=%g, emp margin long put=%g\n", marginLongPut, empMarginLongPut)
			t.Errorf("Error=%g is more than tolerance", errorPut)
		}

	}

}

// TestTimeTakenForOptionsRiskFactor
// At the moment you cannot fail this test.
// Eventually we may impose minimum performance here to check
// that we don't have a regression that really slowed things down
// for one reason or another.
func TestTimeTakenForOptionsRiskFactor(t *testing.T) {
	var paramsBs ModelParamsBS
	paramsBs.mu = 0.1
	paramsBs.r = 0.01
	paramsBs.sigma = 1.2

	const numRuns int = 1000000
	start := time.Now()

	for i := 0; i < numRuns; i++ {
		RiskFactorsCall(0.01, 1.0/365.25/24, 1, 1, 0.5, paramsBs)
		RiskFactorsPut(0.01, 1.0/365.25/24, 1, 1, 0.5, paramsBs)
	}

	elapsed := time.Since(start)

	fmt.Printf("Num of times we can calculate call and put risk factors in BS model: %.1f per second.\n",
		float64(numRuns)/elapsed.Seconds())
}
