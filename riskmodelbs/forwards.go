package riskmodelbs

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

// ModelParamsBS collect the parameters of Black-Scholes model.
// Here mu is the real-world measure growth rate, r is the risk-free interest rate, sigma is volatiliy
type ModelParamsBS struct {
	Mu    float64
	R     float64
	Sigma float64
}

// RiskFactorsForward calculates the risk factors based on Black Scholes model for the evolution
// of the risky asset (i.e. geometric brownian motion i.e. future is lognormal)
func RiskFactorsForward(lambd, tau float64, modelParams ModelParamsBS) RiskFactors {
	mu := modelParams.Mu
	sigma := modelParams.Sigma
	muBar := (mu - 0.5*sigma*sigma) * tau
	sigmaBar := math.Sqrt(tau) * sigma

	riskFactorShort := riskmeasures.NegativeLogNormalEs(muBar, sigmaBar, lambd) - 1.0
	riskFactorLong := riskmeasures.LogNormalEs(muBar, sigmaBar, lambd) + 1.0

	factors := RiskFactors{riskFactorLong, riskFactorShort}
	return factors
}
