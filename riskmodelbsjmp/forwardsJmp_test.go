package riskmodelbsjmp

import (
	"fmt"
	"math"
	"testing"
	"time"

	"code.vegaprotocol.io/quant/riskmodelbs"
)

const testTolerance float64 = 1.0e-8

func TestMartingale(t *testing.T) {
	T := 1.0 / 365.25
	mu := 0.2
	sigma := 2.0
	gamma := 10.0
	jmpMeanA := 0.0
	jmpStddevB := 1.0
	N := 5000000
	tol := 1e-3

	diffPart, jmpPart := generateJumpDiffSamples(T, mu, sigma, gamma, jmpMeanA, jmpStddevB, N)

	// estimate mean of e^{\mu T} \mathbb E[e^{X_T}]
	mean := 0.0
	expMinusMuT := math.Exp(-mu * T)
	for i := 0; i < N; i++ {
		mean += expMinusMuT * diffPart[i] * jmpPart[i]
	}
	mean = mean / float64(N)
	t.Logf("emp mean=%g, should be =1.0\n", mean)
	if math.Abs(1.0-mean) > tol {
		t.Errorf("Failed the martingale test.")
	}
}

// TestBSJmpFwdRiskFactors tests fft method vs. monte carlo
func TestBSJmpFwdRiskFactors(t *testing.T) {
	const tol = 1e-2
	nSamples := 500000 // for MC risk factors

	const lambda float64 = 0.1
	const numRuns int = 1

	for runIdx := 0; runIdx < numRuns; runIdx++ {
		var paramsBsJmp RiskModelParamsBSJmp
		paramsBsJmp.Mu = 0.0
		paramsBsJmp.R = 0.00
		paramsBsJmp.Sigma = 2.0
		paramsBsJmp.Gamma = 10.0
		paramsBsJmp.JmpMeanA = 0.0
		paramsBsJmp.JmpStddevB = 1.0
		tau := 1.0 / 365.25

		//paramsBs.mu = distuv.UnitUniform.Rand()
		//paramsBs.r = distuv.UnitUniform.Rand()
		//paramsBs.sigma = distuv.UnitUniform.Rand()
		//tau := distuv.UnitUniform.Rand()

		riskFactors, err := RiskFactorsForward(lambda, tau, paramsBsJmp)

		if err != nil || math.IsNaN(riskFactors.Short) || math.IsNaN(riskFactors.Long) ||
			math.IsInf(riskFactors.Short, 0) || math.IsInf(riskFactors.Long, 0) {

			t.Logf("r=%g, mu=%g, sigma=%g, tau=%g\n", paramsBsJmp.R, paramsBsJmp.Mu, paramsBsJmp.Sigma, tau)
			t.Errorf("risk factor is NaN or Inf")
		}

		riskFactorsEmpirical, errEmp := EmpiricalRiskFactorsForward(lambda, tau, paramsBsJmp,
			nSamples)
		if errEmp != nil {
			t.Errorf("fail")
			return
		}

		errShort := math.Abs(riskFactors.Short - riskFactorsEmpirical.Short)
		errLong := math.Abs(riskFactors.Long - riskFactorsEmpirical.Long)
		if (errShort > tol) || (errLong > tol) {
			t.Logf("r=%g, mu=%g, sigma=%g, tau=%g\n", paramsBsJmp.R, paramsBsJmp.Mu, paramsBsJmp.Sigma, tau)
			t.Logf("rf short=%g, rf long=%g\n", riskFactors.Short, riskFactors.Long)
			t.Logf("emp rf short=%g, emp rf long=%g\n", riskFactorsEmpirical.Short,
				riskFactorsEmpirical.Long)
			var paramsBs riskmodelbs.RiskModelParamsBS
			paramsBs.Mu = paramsBsJmp.Mu
			paramsBs.Sigma = paramsBsJmp.Sigma
			bsRiskFactors := riskmodelbs.RiskFactorsForward(lambda, tau, paramsBs)
			t.Logf("For comparison only: the corresponding BS risk factors.")
			t.Logf("BS rf short=%g, BS  rf long=%g\n", bsRiskFactors.Short, bsRiskFactors.Long)
			t.Errorf("risk factor mismatch between FFT and MC")
		}
	}
}

func TestTimeTakenForForwardsRiskFactor(t *testing.T) {
	var paramsBsJmp RiskModelParamsBSJmp
	paramsBsJmp.Mu = 0.1
	paramsBsJmp.R = 0.01
	paramsBsJmp.Sigma = 1.2
	paramsBsJmp.Gamma = 1.0
	paramsBsJmp.JmpMeanA = 0.0
	paramsBsJmp.JmpStddevB = 1.0

	const numRuns int = 2
	start := time.Now()

	for i := 0; i < numRuns; i++ {
		RiskFactorsForward(0.01, 1.0/365.25/24, paramsBsJmp)
	}

	elapsed := time.Since(start)

	fmt.Printf("Num of times we can calculate forward risk factors in BS Jmp model: %.1f per second.\n",
		float64(numRuns)/elapsed.Seconds())
}
