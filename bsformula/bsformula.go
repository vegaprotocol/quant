package bsformula

import . "math"
import . "gitlab.com/vega-protocol/quant/misc"

func d1Fn(S, K, r, sigma, T float64) float64 {
	return (Log(S/K) + (r+sigma*sigma*0.5)*T) / (sigma * Sqrt(T))
}

// BSCallProb1 returns the P_1 in call = S P_1 - Ke^(-rT)P_2
func BSCallProb1(S, K, r, sigma, T float64) float64 {
	var d1 = d1Fn(S, K, r, sigma, T)
	return ApproxGaussCdf(d1)
}

// BSCallProb2 returns the P_2 in call = S P_1 - Ke^(-rT)P_2
func BSCallProb2(S, K, r, sigma, T float64) float64 {
	var d1 = d1Fn(S, K, r, sigma, T)
	var d2 = d1 - sigma*Sqrt(T)
	return ApproxGaussCdf(d2)
}

// BSCallPrice calculates the call option price according to the BS formula
func BSCallPrice(S, K, r, sigma, T float64) float64 {
	var d1 = d1Fn(S, K, r, sigma, T)
	var d2 = d1 - sigma*Sqrt(T)
	return S*ApproxGaussCdf(d1) - K*Exp(-r*T)*ApproxGaussCdf(d2)
}

func getPutFromCallPrice(S, K, r, T, callPrice float64) float64 {
	return callPrice - S + K*Exp(-r*T)
}

// BSPutPrice calculates the put option price according to the BS formula
func BSPutPrice(S, K, r, sigma, T float64) float64 {
	return getPutFromCallPrice(S, K, r, T, BSCallPrice(S, K, r, sigma, T))
}

// BSCallDelta calculates the BS Delta (partial derivative w.r.t. S)
func BSCallDelta(S, K, r, sigma, T float64) float64 {
	return ApproxGaussCdf(d1Fn(S, K, r, sigma, T))
}

// BSPutDelta calculates the BS Delta (partial derivative w.r.t. S)
func BSPutDelta(S, K, r, sigma, T float64) float64 {
	return -ApproxGaussCdf(-d1Fn(S, K, r, sigma, T))
}

// BSVega calculates the BS Vega (partial derivative w.r.t. sigma)
// Note that it's identical for both puts and calls
func BSVega(S, K, r, sigma, T float64) float64 {
	d1 := d1Fn(S, K, r, sigma, T)
	return S * GaussDensity(d1) * Sqrt(T)
}

// ImpliedVol calculates the implied volatility
// from call or put price as indicated by isCall
func ImpliedVol(S, K, r, T, price float64, isCall bool) (float64, error) {
	const solverTol float64 = 1e-12
	const solverMaxIt int = 100

	var guessVol float64 = 0.8
	var bsAsFnOfSigma func(sigma float64) float64
	if isCall {
		bsAsFnOfSigma = func(sigma float64) float64 {
			return BSCallPrice(S, K, r, sigma, T) - price
		}
	} else {
		bsAsFnOfSigma = func(sigma float64) float64 {
			return BSPutPrice(S, K, r, sigma, T) - price
		}
	}
	bsAsFnOfSigmaPrime := func(sigma float64) float64 {
		return BSVega(S, K, r, sigma, T)
	}

	impliedVol, err := FindRoot(bsAsFnOfSigma, bsAsFnOfSigmaPrime, guessVol, solverMaxIt, solverTol)
	if err == nil {
		return impliedVol, nil
	}
	return NaN(), err
}
