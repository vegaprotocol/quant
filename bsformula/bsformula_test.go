package bsformula

import (
	"fmt"
	"math"
	"math/rand"
	"testing"
	"time"

	"gonum.org/v1/gonum/stat/distuv"
)

// the inputs are S, K, r, sigma, T
var testValues = []struct {
	S     float64
	K     float64
	r     float64
	sigma float64
	T     float64
}{
	{1.0, 1.0, 0.02, 0.5, 0.25},  // test at the money
	{1.0, 1.0, -0.02, 0.5, 0.25}, // test negative interest rate
	{1.0, 1.5, 0.02, 0.5, 0.25},  // test high strike
	{1.0, 0.5, 0.02, 0.5, 0.25},  // test low strike
	{1.0, 0.5, 0.0, 0.5, 0.25},   // test zero interest rate
	{1.0, 0.5, 0.0, 0.1, 2},      // test long maturity
	{75, 80, 0.0, 0.1, 2},        // test non-unit strike / mat
}

// TestBSPricesUsingMonteCarlo we use MC simulation to check call
// and put prices.
func TestBSPricesUsingMonteCarlo(t *testing.T) {
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

	for _, table := range testValues {
		S := table.S
		K := table.K
		r := table.r
		sigma := table.sigma
		T := table.T

		bsCall := BSCallPrice(S, K, r, sigma, T)
		bsPut := BSPutPrice(S, K, r, sigma, T)
		var callPayoff float64 // = 0.0
		var putPayoff float64  // = 0.0
		for i := 0; i < numMCSamples; i++ {
			SatT := S * math.Exp((r-0.5*sigma*sigma)*T+sigma*math.Sqrt(T)*Z[i])
			callPayoff += math.Max(SatT-K, 0.0)
			putPayoff += math.Max(K-SatT, 0.0)
		}
		mcCall := math.Exp(-r*T) * callPayoff / float64(numMCSamples)
		mcPut := math.Exp(-r*T) * putPayoff / float64(numMCSamples)
		error := (math.Abs(mcCall-bsCall) + math.Abs(mcPut-bsPut)) / S
		if math.IsNaN(error) || math.IsInf(error, 0) {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is NaN or Inf!\n",
				S, K, r, sigma, T, error)
		}

		if error > testToleranceForMC {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is greater than tolerance.\n",
				S, K, r, sigma, T, error)
		}

	}
}

// TestBSProbsVsFormula tests whether bsCallCalcProb1 and bsCallCalcProb2 behave correctly
func TestBSProbsVsFormula(t *testing.T) {

	const testTolerance float64 = 1.0e-16

	for _, table := range testValues {
		S := table.S
		K := table.K
		r := table.r
		sigma := table.sigma
		T := table.T

		bsCall := BSCallPrice(S, K, r, sigma, T)
		bsCallDirect := S*BSCallProb1(S, K, r, sigma, T) - K*math.Exp(-r*T)*BSCallProb2(S, K, r, sigma, T)
		error := math.Abs(bsCall-bsCallDirect) / S

		if math.IsNaN(error) || math.IsInf(error, 0) {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is NaN or Inf!\n",
				S, K, r, sigma, T, error)
		}

		if error > testTolerance {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is greater than tolerance.\n",
				S, K, r, sigma, T, error)
		}

	}
}

// TestBSDeltaVsFiniteDifference we use finite difference to check the BS deltas
func TestBSDeltaVsFiniteDifference(t *testing.T) {
	const testTolerance float64 = 1.0e-3

	const delta float64 = 1.0e-9

	for _, table := range testValues {
		S := table.S
		K := table.K
		r := table.r
		sigma := table.sigma
		T := table.T

		callDelta := BSCallDelta(S, K, r, sigma, T)
		putDelta := BSPutDelta(S, K, r, sigma, T)

		bsCallPlusBump := BSCallPrice(S+delta, K, r, sigma, T)
		bsPutPlusBump := BSPutPrice(S+delta, K, r, sigma, T)

		bsCallMinusBump := BSCallPrice(S-delta, K, r, sigma, T)
		bsPutMinusBump := BSPutPrice(S-delta, K, r, sigma, T)

		callDeltaFd := (bsCallPlusBump - bsCallMinusBump) / (2 * delta)
		putDeltaFd := (bsPutPlusBump - bsPutMinusBump) / (2 * delta)

		error := (math.Abs(callDelta-callDeltaFd) + math.Abs(putDelta-putDeltaFd)) / S

		if math.IsNaN(error) || math.IsInf(error, 0) {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is NaN or Inf!\n",
				S, K, r, sigma, T, error)
		}

		if error > testTolerance {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is greater than tolerance.\n",
				S, K, r, sigma, T, error)
		}
	}
}

func TestBSVegaVsFiniteDifference(t *testing.T) {
	const testTolerance float64 = 1.0e-3
	const delta float64 = 1.0e-9

	for _, table := range testValues {
		S := table.S
		K := table.K
		r := table.r
		sigma := table.sigma
		T := table.T

		vega := BSVega(S, K, r, sigma, T)

		bsCallPlusBump := BSCallPrice(S, K, r, sigma+delta, T)
		bsPutPlusBump := BSPutPrice(S, K, r, sigma+delta, T)

		bsCallMinusBump := BSCallPrice(S, K, r, sigma-delta, T)
		bsPutMinusBump := BSPutPrice(S, K, r, sigma-delta, T)

		callVegaFd := (bsCallPlusBump - bsCallMinusBump) / (2 * delta)
		putVegaFd := (bsPutPlusBump - bsPutMinusBump) / (2 * delta)

		error := (math.Abs(vega-callVegaFd) + math.Abs(vega-putVegaFd)) / S

		if math.IsNaN(error) || math.IsInf(error, 0) {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is NaN or Inf!\n",
				S, K, r, sigma, T, error)
		}

		if error > testTolerance {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is greater than tolerance.\n",
				S, K, r, sigma, T, error)
		}
	}
}

func TestImpliedVolCalcs(t *testing.T) {
	const testTolerance float64 = 1.0e-5

	// the inputs are S, K, r, sigma, T
	testValuesImpVol := []struct {
		S     float64
		K     float64
		r     float64
		sigma float64
		T     float64
	}{
		{1.0, 1.0, 0.02, 0.5, 0.25},  // test at the money
		{1.0, 1.0, -0.02, 0.5, 0.25}, // test negative interest rate
		{1.0, 1.5, 0.02, 0.5, 0.25},  // test high strike
		{1.0, 0.5, 0.02, 0.5, 0.25},  // test low strike
		{1.0, 0.5, 0.0, 0.5, 0.25},   // test zero interest rate
		{1.0, 0.5, 0.0, 0.1, 2},      // test long maturity
		{75, 80, 0.0, 0.1, 2},        // test non-unit strike / mat
	}

	for _, table := range testValuesImpVol {
		S := table.S
		K := table.K
		r := table.r
		sigma := table.sigma
		T := table.T

		call := BSCallPrice(S, K, r, sigma, T)
		put := BSPutPrice(S, K, r, sigma, T)

		volFromCall, errCall := ImpliedVol(S, K, r, T, call, true)
		if errCall != nil {
			t.Errorf(errCall.Error())
		}

		volFromPut, errPut := ImpliedVol(S, K, r, T, put, false)
		if errPut != nil {
			t.Errorf(errPut.Error())
		}

		error := (math.Abs(volFromCall-sigma) + math.Abs(volFromPut-sigma))

		if math.IsNaN(error) || math.IsInf(error, 0) {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is NaN or Inf!\n",
				S, K, r, sigma, T, error)
		}

		if error > testTolerance {
			t.Errorf("S=%g, K=%g, r=%g, sigma=%g, T=%g, error=%g is greater than tolerance.\n",
				S, K, r, sigma, T, error)
		}
	}

	// Test what happens if solver fails
	impVol, err := ImpliedVol(100, 100, 0, 1, 0, true)
	if err == nil {
		t.Errorf("This should have resulted in error.")
		t.Logf("Claimed implied vol=%g", impVol)
	}
}

func TestTimeTakenForBSFormula(t *testing.T) {
	const S float64 = 1
	const K float64 = 1
	const sigma float64 = 1.5
	const mu float64 = 0.1
	const tau float64 = 1.0 / 365.25 / 24

	const numMCSamples int = 10000
	const timingN int = 1000
	start := time.Now()
	for timingIdx := 0; timingIdx < timingN; timingIdx++ {
		callPriceAfterTau := make([]float64, numMCSamples)
		for i := 0; i < numMCSamples; i++ {
			z := rand.NormFloat64()
			Stau := S * math.Exp((mu-0.5*sigma*sigma)*tau+math.Sqrt(tau)*sigma*z)
			callPriceAfterTau[i] = BSCallPrice(Stau, K, 0.05, sigma, 1.0)
		}
	}
	elapsed := time.Since(start)

	fmt.Printf("Num of times we can generate %v samples from GBM and calc. BS price for each: %.1f per second.\n",
		numMCSamples, float64(timingN)/elapsed.Seconds())
}
