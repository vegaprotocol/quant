package riskmodelbsjmp

import (
	"math"

	"code.vegaprotocol.io/quant/riskmeasures"
)

// RiskFactors are used in margin calculation as follows:
// the margin is set to unrealised P&L + mark price x risk factor x open volume;
// use the Long one if the overall position of the participant is long (i.e open volume > 0)
// and use Short if open volume < 0
type RiskFactors struct {
	Long  float64
	Short float64
}

// RiskModelParamsBSJmp collect the parameters of Black-Scholes model.
// Here mu is the real-world measure growth rate, r is the risk-free interest rate, sigma is volatiliy
// gamma is jump intensity (jump inter-arrival times are exponential with this param)
// finally the jump sizes are normally distributed N(a,b^2),
// jmpMeanA = a, variance = jmpStddevB^2 b^2
type RiskModelParamsBSJmp struct {
	Mu         float64 // real world growth
	R          float64 // interest rate
	Sigma      float64 // volatility of diffusion part
	Gamma      float64 // jump intensity (jump inter-arrival times are exponential with this param)
	JmpMeanA   float64 // jump mean
	JmpStddevB float64 // jump stddev
}

// RiskFactorsForward calculates the risk factors based on Black Scholes model for the evolution
// of the risky asset (i.e. geometric brownian motion i.e. future is lognormal)
func RiskFactorsForward(lambd, tau float64, modelParams RiskModelParamsBSJmp) (RiskFactors, error) {

	negativeLogNormalWithJmpEs, err1 := riskmeasures.NegativeLogNormalJmpEs(tau,
		modelParams.Mu,
		modelParams.Sigma,
		modelParams.Gamma,
		modelParams.JmpMeanA,
		modelParams.JmpStddevB,
		lambd)
	if err1 != nil {
		factors := RiskFactors{math.NaN(), math.NaN()}
		return factors, err1
	}

	logNormalWithJmpEs, err2 := riskmeasures.LogNormalJmpEs(tau,
		modelParams.Mu,
		modelParams.Sigma,
		modelParams.Gamma,
		modelParams.JmpMeanA,
		modelParams.JmpStddevB,
		lambd)
	if err2 != nil {
		factors := RiskFactors{math.NaN(), math.NaN()}
		return factors, err2
	}

	riskFactorShort := negativeLogNormalWithJmpEs - 1.0
	riskFactorLong := logNormalWithJmpEs + 1.0

	factors := RiskFactors{riskFactorLong, riskFactorShort}
	return factors, nil
}

// EmpiricalRiskFactorsForward ...
func EmpiricalRiskFactorsForward(lambd, tau float64, modelParams RiskModelParamsBSJmp, numSamples int) (RiskFactors, error) {
	expJumpDiffSamples := make([]float64, numSamples)
	minusExpJumpDiffSamples := make([]float64, numSamples)
	for i := 0; i < numSamples; i++ {
		diff, jmp := generateJumpDiffSample(tau, modelParams.Mu, modelParams.Sigma,
			modelParams.Gamma, modelParams.JmpMeanA, modelParams.JmpStddevB)
		expJumpDiffSamples[i] = math.Exp(diff + jmp)
		minusExpJumpDiffSamples[i] = -expJumpDiffSamples[i]
	}
	logNormalWithJmpEs := riskmeasures.EmpiricalEs(expJumpDiffSamples, lambd, false)
	negativeLogNormalWithJmpEs := riskmeasures.EmpiricalEs(minusExpJumpDiffSamples, lambd, false)
	riskFactorShort := negativeLogNormalWithJmpEs - 1.0
	riskFactorLong := logNormalWithJmpEs + 1.0

	factors := RiskFactors{riskFactorLong, riskFactorShort}
	return factors, nil
}
