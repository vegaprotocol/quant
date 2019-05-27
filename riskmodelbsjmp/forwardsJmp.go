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
	mu         float64 // real world growth
	r          float64 // interest rate
	sigma      float64 // volatility of diffusion part
	gamma      float64 // jump intensity (jump inter-arrival times are exponential with this param)
	jmpMeanA   float64 // jump mean
	jmpStddevB float64 // jump stddev
}

// RiskFactorsForward calculates the risk factors based on Black Scholes model for the evolution
// of the risky asset (i.e. geometric brownian motion i.e. future is lognormal)
func RiskFactorsForward(lambd, tau float64, modelParams RiskModelParamsBSJmp) (RiskFactors, error) {

	negativeLogNormalWithJmpEs, err1 := riskmeasures.NegativeLogNormalJmpEs(tau,
		modelParams.mu,
		modelParams.sigma,
		modelParams.gamma,
		modelParams.jmpMeanA,
		modelParams.jmpStddevB,
		lambd)
	if err1 != nil {
		factors := RiskFactors{math.NaN(), math.NaN()}
		return factors, err1
	}

	logNormalWithJmpEs, err2 := riskmeasures.LogNormalJmpEs(tau,
		modelParams.mu,
		modelParams.sigma,
		modelParams.gamma,
		modelParams.jmpMeanA,
		modelParams.jmpStddevB,
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
