package riskmeasures

import (
	"math"

	"code.vegaprotocol.io/quant/fftpricing"
)

// LogNormalJmpEs returns the expected shortfall of exponential of jump diffusion.
// The inputs are the model parameters; the method uses the expression of expected
// shortfall in terms of minimization of expectation of a call-option-like expression;
// the call option is valued using FFT method.
func LogNormalJmpEs(T, mu, sigma, gamma, a, b, lambd float64) (float64, error) {
	// set up parameters
	var p fftpricing.PricingModelParamsBSJmp
	p.R = mu
	p.Sigma = sigma
	p.Gamma = gamma
	p.JmpMeanA = a
	p.JmpStddevB = b

	// we'll use later
	expMuT := math.Exp(mu * T)
	// get call prices for all strikes
	ks, prices, err := fftpricing.JumpDiffCallPricesFftFullAxis(T, p)
	if err != nil {
		return math.NaN(), err
	}
	minVal := 0.0
	for i := 0; i < len(ks); i++ {
		k := ks[i]
		expK := math.Exp(k)
		newVal := -expK + (1.0/lambd)*(expK-expMuT+expMuT*prices[i])
		minVal = math.Min(minVal, newVal)
	}
	return -minVal, nil
}

// NegativeLogNormalJmpEs returns the expected shortfall of
// minus one times exponential of jump diffusion.
// The inputs are the model parameters; the method uses the expression of expected
// shortfall in terms of minimization of expectation of a call-option-like expression;
// the call option is valued using FFT method.
func NegativeLogNormalJmpEs(T, mu, sigma, gamma, a, b, lambd float64) (float64, error) {
	// set up parameters
	var p fftpricing.PricingModelParamsBSJmp
	p.R = mu
	p.Sigma = sigma
	p.Gamma = gamma
	p.JmpMeanA = a
	p.JmpStddevB = b

	// we'll use later
	expMuT := math.Exp(mu * T)
	// get call prices for all strikes
	ks, prices, err := fftpricing.JumpDiffCallPricesFftFullAxis(T, p)
	if err != nil {
		return math.NaN(), err
	}
	minVal := (1.0 / lambd) * expMuT
	for i := 0; i < len(ks); i++ {
		k := ks[i]
		expK := math.Exp(k)
		newVal := expK + (1.0/lambd)*expMuT*prices[i]
		minVal = math.Min(minVal, newVal)
	}
	return -minVal, nil
}
