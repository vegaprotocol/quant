package riskmeasures

import (
	"math"
	"sort"
	"sync"
	"testing"

	"golang.org/x/exp/rand"
)

const testTolerance float64 = 1.0e-8

func TestVarAndESLogNormalAgainsPrecomputed(t *testing.T) {
	tables := []struct {
		mu     float64
		sigma  float64
		lambda float64
		VaR    float64
		ES     float64
	}{
		{0.0, 0.1, 0.01, -7.9244293082e-01, -7.6640460826e-01},
		{0.0, 0.1, 0.05, -8.4833017412e-01, -8.1416417294e-01},
		{0.0, 1.0, 0.01, -9.7651733070e-02, -7.2537170781e-02},
		{1.0, 2.0, 0.01, -2.5921157598e-02, -1.5225031293e-02},
	}

	for _, table := range tables {
		VaR := LogNormalVaR(table.mu, table.sigma, table.lambda)
		errorVaR := math.Abs(VaR - table.VaR)
		if math.IsNaN(errorVaR) || math.IsInf(errorVaR, 0) || errorVaR > testTolerance {
			t.Errorf("Test VaR=%g, precomputed VaR=%g, Error=%g is greater than tolerance.\n", errorVaR, table.VaR, errorVaR)
		}
		ES := LogNormalEs(table.mu, table.sigma, table.lambda)
		errorES := math.Abs(ES - table.ES)
		if math.IsNaN(errorES) || math.IsInf(errorES, 0) || errorES > testTolerance {
			t.Errorf("Test ES=%g, precomputed ES=%g, Error=%g is greater than tolerance.\n", errorES, table.ES, errorES)
		}
	}
}

// Table of values which will be used for
// testing.
var testValsForESLognormal = []struct {
	mu     float64
	sigma  float64
	lambda float64
}{
	{0.0, 0.1, 0.01},
	{0.0, 0.1, 0.05},
	{0.0, 0.5, 0.01},
	{1.0, 0.6, 0.03},
	{-1.0, 0.3, 0.01},
	{-1.0, 1.2, 0.2},
	{-1.0, 0.8, 0.5},
	{-1.0, 0.5, 0.7},
}

func runTestVaRAndESLogNormalUsingMC(t *testing.T, positiveLogNormal bool) {
	// these two constants are related, if the test takes too long then reduce num samples but
	// increase tolerance... (but this decreases the strength of the test I am afraid)
	const testToleranceForMC float64 = 2.5e-3
	const numMCSamples int = 8000000

	Z := make([]float64, 2*numMCSamples)
	for i := 0; i < numMCSamples; i++ {
		//Z[i] = distuv.UnitNormal.Rand()
		z := rand.NormFloat64()
		Z[i] = z
		Z[numMCSamples+i] = -z //antithetic
	}
	// since the exponential transformation is monotone and sigma > 0
	// we can assume that once Z has been sorted exp(mu + sigma*Z) is also sorted
	sort.Float64s(Z)
	var wg sync.WaitGroup
	for _, table := range testValsForESLognormal {
		mu := table.mu
		sigma := table.sigma
		lambda := table.lambda

		wg.Add(1)
		go func(mu, sigma, lambda float64) {
			defer wg.Done()
			var varExact float64
			var sign float64
			var esExact float64
			if positiveLogNormal {
				varExact = LogNormalVaR(mu, sigma, lambda)
				esExact = LogNormalEs(mu, sigma, lambda)
				sign = 1.0
			} else {
				varExact = NegativeLogNormalVaR(mu, sigma, lambda)
				esExact = NegativeLogNormalEs(mu, sigma, lambda)
				sign = -1.0
			}

			// generate lognormal samples
			X := make([]float64, 2*numMCSamples)
			for i := 0; i < 2*numMCSamples; i++ {
				X[i] = sign * math.Exp(mu+sigma*Z[i])
			}
			varEmpirical := EmpiricalVaR(X, lambda, false) // to cover the test case
			errorVaR := math.Abs(varEmpirical - varExact)
			if math.IsNaN(errorVaR) || math.IsInf(errorVaR, 0) || errorVaR > testToleranceForMC {
				t.Errorf("VaR: Error greater than MC tolerance (%g): mu=%g, sigma=%g, lambda=%g, exact=%g, empirical=%g, error=%g\n", testToleranceForMC, mu, sigma, lambda, varExact, varEmpirical, errorVaR)
			}

			esEmpirical := EmpiricalEs(X, lambda, true)
			errorES := math.Abs(esExact - esEmpirical)
			if math.IsNaN(errorES) || math.IsInf(errorES, 0) || errorES > testToleranceForMC {
				t.Errorf("ES : Error greater than MC tolerance (%g): mu=%g, sigma=%g, lambda=%g, exact=%g, empirical=%g, error=%g\n", testToleranceForMC, mu, sigma, lambda, varExact, varEmpirical, errorES)
			}
		}(mu, sigma, lambda)
	}
	wg.Wait()
}

func TestVaRAndESLogNormalUsingMC(t *testing.T) {
	runTestVaRAndESLogNormalUsingMC(t, true)
	runTestVaRAndESLogNormalUsingMC(t, false)
}

func TestVarAndESNegLogNormalAgainsPrecomputed(t *testing.T) {

	tables := []struct {
		mu     float64
		sigma  float64
		lambda float64
		VaR    float64
		ES     float64
	}{
		{0.0, 0.1, 0.01, 1.2619205259e+00, 1.3060584483e+00},
		{0.0, 0.1, 0.05, 1.1787863152e+00, 1.2299511311e+00},
		{0.0, 1.0, 0.01, 1.0240473656e+01, 1.5227960301e+01},
		{1.0, 2.0, 0.01, 2.8505887791e+02, 7.4734383372e+02},
	}

	for _, table := range tables {
		VaR := NegativeLogNormalVaR(table.mu, table.sigma, table.lambda)
		errorVaR := math.Abs(VaR - table.VaR)
		if math.IsNaN(errorVaR) || math.IsInf(errorVaR, 0) || errorVaR > testTolerance {
			t.Errorf("Test VaR=%g, precomputed VaR=%g, Error=%g is greater than tolerance.\n", errorVaR, table.VaR, errorVaR)
		}
		ES := NegativeLogNormalEs(table.mu, table.sigma, table.lambda)
		errorES := math.Abs(ES - table.ES)
		if errorES > testTolerance {
			t.Errorf("Test ES=%g, precomputed ES=%g, Error=%g is greater than tolerance.\n", errorES, table.ES, errorES)
		}
	}
}
