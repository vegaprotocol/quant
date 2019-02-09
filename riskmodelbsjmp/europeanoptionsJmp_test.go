package riskmodelbsjmp

import (
	"testing"
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

}

// TestTimeTakenForOptionsRiskFactor
// At the moment you cannot fail this test.
// Eventually we may impose minimum performance here to check
// that we don't have a regression that really slowed things down
// for one reason or another.
func TestTimeTakenForOptionsRiskFactor(t *testing.T) {

}
