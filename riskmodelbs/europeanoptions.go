package riskmodelbs

import (
	"math"

	"code.vegaprotocol.io/quant/bsformula"
	"code.vegaprotocol.io/quant/riskmeasures"
)

// RiskFactorsCall calculates the risk factors based on Black Scholes model for the evolution
// of the risky asset (i.e. geometric brownian motion i.e. risky asset dist. is lognormal)
// The risk factors returned are for CALL option
func RiskFactorsCall(lambd, tau, S, K, T float64, p ModelParamsBS) RiskFactors {
	muBar := (p.mu - 0.5*p.sigma*p.sigma) * tau
	sigmaBar := math.Sqrt(tau) * p.sigma

	bsProb1 := bsformula.BSCallProb1(S, K, p.r, p.sigma, T)
	negLogNormEs := riskmeasures.NegativeLogNormalEs(muBar, sigmaBar, lambd)
	riskFactorShort := bsProb1 * (negLogNormEs - 1.0)

	logNormEs := riskmeasures.LogNormalEs(muBar, sigmaBar, lambd)
	riskFactorLong := bsProb1 * (logNormEs + 1.0)

	factors := RiskFactors{riskFactorLong, riskFactorShort}
	return factors
}

// RiskFactorsPut calculates the risk factors based on Black Scholes model for the evolution
// of the risky asset (i.e. geometric brownian motion i.e. risky asset dist. is lognormal)
// The risk factors returned are for PUT option
func RiskFactorsPut(lambd, tau, S, K, T float64, p ModelParamsBS) RiskFactors {
	muBar := (p.mu - 0.5*p.sigma*p.sigma) * tau
	sigmaBar := math.Sqrt(tau) * p.sigma

	bsProb1 := bsformula.BSCallProb1(S, K, p.r, p.sigma, T)

	logNormEs := riskmeasures.LogNormalEs(muBar, sigmaBar, lambd)
	riskFactorShort := (1.0 - bsProb1) * (logNormEs + 1.0)

	negLogNormEs := riskmeasures.NegativeLogNormalEs(muBar, sigmaBar, lambd)
	riskFactorLong := (1.0 - bsProb1) * (negLogNormEs - 1.0)

	factors := RiskFactors{riskFactorLong, riskFactorShort}
	return factors
}
