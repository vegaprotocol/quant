package riskmodelbsjmp

import (
	"math"
	"testing"

	"code.vegaprotocol.io/quant/fftpricing"
)

func TestFFTFullAxisMethodForJmpAllStrikes(t *testing.T) {
	const numMCsamples int = 20000
	const tol float64 = 1e-2
	const minStrike float64 = 0.5 // we ignore small strikes
	const maxStrike float64 = 1.5 // we ignore large strikes

	const T float64 = 1.0 / 365.25

	paramsJmp := fftpricing.CreatePricingParamsForBSJmp(0.0, 2., 100., 0., 0.0000001)
	//paramsJmp.r = 0.0
	//paramsJmp.sigma = 0.5
	//paramsJmp.gamma = 0.1
	//paramsJmp.jmpMeanA = 0.4
	//paramsJmp.jmpStddevB = 0.05

	strikes, optValsFft, err := fftpricing.JumpDiffCallPricesFftFullAxis(T, paramsJmp)
	if err != nil {
		t.Errorf(err.Error())
		return
	}

	N := len(optValsFft)

	for i := 0; i < N; i++ {
		strike := strikes[i]
		// we ignore small and large strikes
		if strike < maxStrike && strike > minStrike {
			mcPrice := EuropeanCallOptionPriceMC(1.0, strike, T, paramsJmp, numMCsamples)
			fftPrice := optValsFft[i]
			err := math.Abs(mcPrice - fftPrice)

			if err > tol && strike < 1.5 && strike > 0.5 {
				t.Errorf("err = %g, MC = %g, FFT = %g at K = %g", err, mcPrice, fftPrice, strikes[i])
				return
			}
		}
	}
}
