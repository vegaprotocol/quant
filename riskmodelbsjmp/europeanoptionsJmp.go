package riskmodelbsjmp

import (
	"math"

	"code.vegaprotocol.io/quant/bsformula"
	"code.vegaprotocol.io/quant/riskmeasures"
	"gonum.org/v1/gonum/stat/distuv"
)

func generateJumpDiffSample(T, mu, sigma, gamma, jmpMeanA, jmpStddevB float64) float64 {
	// I am not sure whether I need to set the random source here as well
	var poissonParams distuv.Poisson
	poissonParams.Lambda = gamma * T

	// normal part
	z1 := distuv.UnitNormal.Rand()
	x1 := (mu-0.5*sigma*sigma)*T + math.Sqrt(T)*sigma*z1

	// jump part
	z2 := distuv.UnitNormal.Rand()
	alpha := math.Exp(jmpMeanA+0.5*jmpStddevB*jmpStddevB) - 1.0

	// how many jumps
	Njumps := distuv.Poisson.Rand(poissonParams)
	x2 := jmpMeanA*Njumps + jmpStddevB*math.Sqrt(Njumps)*z2 - gamma*alpha*T
	XT := x1 + x2

	return XT
}

// EuropeanCallOptionPriceMC uses direct MC simulation to approximate the
// call price in the jump diffusion model
func EuropeanCallOptionPriceMC(S, K, T float64, p ModelParamsBSJmp) float64 {
	const N int = 10000

	runningSum := 0.0
	for i := 0; i < N; i++ {
		ST := math.Exp(generateJumpDiffSample(T, p.r, p.sigma, p.gamma, p.jmpMeanA, p.jmpStddevB))
		callPayoff := math.Max(ST-K, 0)
		runningSum += callPayoff
	}
	priceMC := math.Exp(-p.r*T) * runningSum / float64(N)
	return priceMC
}

// RiskFactorsCall calculates the risk factors based on Black Scholes model for the evolution
// of the risky asset (i.e. geometric brownian motion i.e. risky asset dist. is lognormal)
// The risk factors returned are for CALL option
func RiskFactorsCall(lambd, tau, S, K, T float64, p ModelParamsBSJmp) RiskFactors {
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
func RiskFactorsPut(lambd, tau, S, K, T float64, p ModelParamsBSJmp) RiskFactors {
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
