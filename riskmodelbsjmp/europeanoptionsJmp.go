package riskmodelbsjmp

import (
	"math"

	"code.vegaprotocol.io/quant/bsformula"
	"code.vegaprotocol.io/quant/fftpricing"
	"code.vegaprotocol.io/quant/misc"
	"code.vegaprotocol.io/quant/riskmeasures"
	"gonum.org/v1/gonum/stat/distuv"
)

// returns the diffusion part and jump part separately
func generateJumpDiffSample(T, mu, sigma, gamma, jmpMeanA, jmpStddevB float64) (float64, float64) {

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

	return x1, x2
}

// returns the diffusion part and jump part separately
func generateJumpDiffSamples(T, mu, sigma, gamma, jmpMeanA, jmpStddevB float64, N int) ([]float64, []float64) {
	samplesSTDiffPt := make([]float64, N)
	samplesSTJumpPt := make([]float64, N)
	for i := 0; i < N; i++ {
		diffPart, jumpPart := generateJumpDiffSample(T, mu, sigma, gamma, jmpMeanA, jmpStddevB)
		samplesSTDiffPt[i] = math.Exp(diffPart)
		samplesSTJumpPt[i] = math.Exp(jumpPart)
	}
	return samplesSTDiffPt, samplesSTJumpPt
}

// EuropeanCallOptionPriceMC uses direct MC simulation to approximate the
// call price in the jump diffusion model
func EuropeanCallOptionPriceMC(S, K, T float64, p fftpricing.PricingModelParamsBSJmp, N int) float64 {

	callPayoffsJmp := make([]float64, N)
	callPayoffsDiff := make([]float64, N)
	for i := 0; i < N; i++ {
		diffPart, jumpPart := generateJumpDiffSample(T, p.R, p.Sigma, p.Gamma, p.JmpMeanA, p.JmpStddevB)
		STJmp := S * math.Exp(diffPart+jumpPart)
		callPayoffsJmp[i] = math.Exp(-p.R*T) * math.Max(STJmp-K, 0)
		STDiff := S * math.Exp(diffPart)
		callPayoffsDiff[i] = math.Exp(-p.R*T) * math.Max(STDiff-K, 0)
	}
	bsPriceAsControl := bsformula.BSCallPrice(1, K, p.R, p.Sigma, T)
	return misc.CalculateControlVariateEstimator(callPayoffsJmp, callPayoffsDiff,
		bsPriceAsControl)
}

// EuropeanCallOptionMCAllStrikes uses direct MC simulation to approximate the
// call price in the jump diffusion model for S=1 and all input strikes
func EuropeanCallOptionMCAllStrikes(strikes []float64, T float64, p fftpricing.PricingModelParamsBSJmp, N int) []float64 {

	samplesSTDiff, samplesSTJump := generateJumpDiffSamples(T, p.R, p.Sigma, p.Gamma, p.JmpMeanA, p.JmpStddevB, N)

	pureDiffPayoffs := make([]float64, N)
	jumpDiffPayoffs := make([]float64, N)

	numStrikes := len(strikes)
	prices := make([]float64, numStrikes)
	for i := 0; i < numStrikes; i++ {
		K := strikes[i]
		for j := 0; j < N; j++ {
			pureDiffPayoffs[j] = math.Exp(-p.R*T) * math.Max(samplesSTDiff[j]-K, 0)
			jumpDiffPayoffs[j] = math.Exp(-p.R*T) * math.Max(samplesSTDiff[j]*samplesSTJump[j]-K, 0)
		}
		bsPriceAsControl := bsformula.BSCallPrice(1, K, p.R, p.Sigma, T)
		prices[i] = misc.CalculateControlVariateEstimator(jumpDiffPayoffs, pureDiffPayoffs, bsPriceAsControl)
	}

	return prices
}

// RiskFactorsCall calculates the risk factors based on Black Scholes model for the evolution
// of the risky asset (i.e. geometric brownian motion i.e. risky asset dist. is lognormal)
// The risk factors returned are for CALL option
func RiskFactorsCall(lambd, tau, S, K, T float64, p RiskModelParamsBSJmp) RiskFactors {
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
func RiskFactorsPut(lambd, tau, S, K, T float64, p RiskModelParamsBSJmp) RiskFactors {
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
