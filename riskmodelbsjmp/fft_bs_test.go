package riskmodelbsjmp

import (
	"fmt"
	"math"
	"math/cmplx"
	"math/rand"
	"testing"
	"time"

	"code.vegaprotocol.io/quant/bsformula"
	"code.vegaprotocol.io/quant/misc"
)

// Fourier transform method for BlackScholes - for testing purposes only price S = 1 and "all" strikes
func bsCallPriceFftHalfAxis(r, sigma, T float64) ([]float64, []float64, error) {

	halfSigmaSquareT := complex(0.5*sigma*sigma*T, 0)
	// characteristic of normal r.v. with variance given by sigma^2 * T
	phiTwo := func(u complex128) complex128 {
		return cmplx.Exp(-halfSigmaSquareT*u*i1 - halfSigmaSquareT*u*u)
	}

	// this is the relevant characteristic function
	zeta := func(v float64) complex128 {
		vC := complex(v, 0)
		return cmplx.Exp(i1*complex(v*r*T, 0)) * (phiTwo(vC-i1) - 1.0) / (i1 * vC * (1 + i1*vC))
	}
	return CallPriceFftHalfAxis(zeta, r, T)
}

// Fourier transform price S = 1 and "all" strikes
func bsCallPriceFftFullAxis(r, sigma, T float64) ([]float64, []float64, error) {

	halfSigmaSquareT := complex(0.5*sigma*sigma*T, 0)
	// characteristic of normal r.v. with variance given by sigma^2 * T
	phiTwo := func(u complex128) complex128 {
		return cmplx.Exp(-halfSigmaSquareT*u*i1 - halfSigmaSquareT*u*u)
	}

	// this is the relevant characteristic function
	zeta := func(v float64) complex128 {
		vC := complex(v, 0)
		return cmplx.Exp(i1*complex(v*r*T, 0)) * (phiTwo(vC-i1) - 1.0) / (i1 * vC * (1 + i1*vC))
	}
	return CallPriceFftFullAxis(zeta, r, T)
}

// TestFFTMethodForBlackScholes can calculate the price of call option in the Black-Scholes
// model using FFT; so if the FFT code is correct then
// the price should closely match the BS formula
func TestFFTMethodForBlackScholesOneValue(t *testing.T) {
	const tol float64 = 1e-3
	const S float64 = 1

	const sigma float64 = 0.5
	const r float64 = 0.05
	const T float64 = 1

	strikesGrid1, optVals1, err1 := bsCallPriceFftHalfAxis(r, sigma, T)
	if err1 != nil {
		t.Errorf(err1.Error())
		return
	}

	bsValsFFT1 := misc.Function{
		X: strikesGrid1,
		Y: optVals1,
	}
	K := 0.8
	bsAtmPrice := bsformula.BSCallPrice(S, K, r, sigma, T)
	fftAtmPrice1 := bsValsFFT1.At(K)
	if math.Abs(bsAtmPrice-fftAtmPrice1) > tol {
		t.Errorf("BS = %g, FFT half axis method = %g", bsAtmPrice, fftAtmPrice1)
	}

	strikesGrid2, optVals2, err2 := bsCallPriceFftFullAxis(r, sigma, T)
	if err2 != nil {
		t.Errorf(err2.Error())
		return
	}

	bsValsFFT2 := misc.Function{
		X: strikesGrid2,
		Y: optVals2,
	}
	fftAtmPrice2 := bsValsFFT2.At(K)
	if math.Abs(bsAtmPrice-fftAtmPrice2) > tol {
		t.Errorf("BS = %g, FFT full axis method = %g", bsAtmPrice, fftAtmPrice2)
	}

}

func TestFFTHalfAxisMethodForBlackScholesAllStrikes(t *testing.T) {
	// the very out of money and very in-the-money values are not that accurate
	const tol float64 = 4e-3
	const sigma float64 = 0.5
	const r float64 = 0.0
	const T float64 = 1

	strikes, optVals, err := bsCallPriceFftHalfAxis(r, sigma, T)
	if err != nil {
		t.Errorf(err.Error())
		return
	}
	N := len(optVals)

	for i := N - 1; i > 0; i-- {
		bsPrice := bsformula.BSCallPrice(1.0, strikes[i], r, sigma, T)
		fftPrice := optVals[i]
		err := math.Abs(bsPrice - fftPrice)
		if err > tol {
			t.Errorf("err = %g, BS = %g, FFT = %g at K = %g", err, bsPrice, fftPrice, strikes[i])
			return
		}
	}
}

func TestFFTFullAxisMethodForBlackScholesAllStrikes(t *testing.T) {
	// the very out of money and very in-the-money values are not that accurate
	const tol float64 = 2.8e-5
	const sigma float64 = 0.5
	const r float64 = 0.0
	const T float64 = 1

	strikes, optVals, err := bsCallPriceFftFullAxis(r, sigma, T)
	if err != nil {
		t.Errorf(err.Error())
		return
	}
	N := len(optVals)

	for i := N - 1; i > 0; i-- {
		bsPrice := bsformula.BSCallPrice(1.0, strikes[i], r, sigma, T)
		fftPrice := optVals[i]
		err := math.Abs(bsPrice - fftPrice)
		if err > tol {
			t.Errorf("err = %g, BS = %g, FFT = %g at K = %g", err, bsPrice, fftPrice, strikes[i])
			return
		}
	}
}

func TestTimeTakenForFFTHalfAxis(t *testing.T) {
	const S float64 = 1
	const K float64 = 1
	const sigma float64 = 0.5
	const T float64 = 1.0
	const mu float64 = 0.1
	const tau float64 = 1.0 / 365.25 / 24

	const numMCSamples int = 10000
	const timingN int = 50
	start := time.Now()
	for timingIdx := 0; timingIdx < timingN; timingIdx++ {
		callPriceAfterTau := make([]float64, numMCSamples)

		strikesGrid, optVals, err := bsCallPriceFftHalfAxis(0.0, sigma, T)
		if err != nil {
			t.Errorf(err.Error())
			return
		}
		bsValsFFT := misc.Function{
			X: strikesGrid,
			Y: optVals,
		}

		for i := 0; i < numMCSamples; i++ {
			z := rand.NormFloat64()
			Stau := S * math.Exp((mu-0.5*sigma*sigma)*tau+math.Sqrt(tau)*sigma*z)
			callPriceAfterTau[i] = bsValsFFT.At(K / Stau)
		}
	}
	elapsed := time.Since(start)

	fmt.Printf("Num of times we can generate %v samples from GBM and calc. opt price via FFT for each: %.1f per second.\n",
		numMCSamples, float64(timingN)/elapsed.Seconds())
}

func TestTimeTakenForFFTFullAxis(t *testing.T) {
	const S float64 = 1
	const K float64 = 1
	const sigma float64 = 0.5
	const T float64 = 1.0
	const mu float64 = 0.1
	const tau float64 = 1.0 / 365.25 / 24

	const numMCSamples int = 10000
	const timingN int = 50
	start := time.Now()
	for timingIdx := 0; timingIdx < timingN; timingIdx++ {
		callPriceAfterTau := make([]float64, numMCSamples)

		strikesGrid, optVals, err := bsCallPriceFftFullAxis(0.0, sigma, T)
		if err != nil {
			t.Errorf(err.Error())
			return
		}
		bsValsFFT := misc.Function{
			X: strikesGrid,
			Y: optVals,
		}

		for i := 0; i < numMCSamples; i++ {
			z := rand.NormFloat64()
			Stau := S * math.Exp((mu-0.5*sigma*sigma)*tau+math.Sqrt(tau)*sigma*z)
			callPriceAfterTau[i] = bsValsFFT.At(K / Stau)
		}
	}
	elapsed := time.Since(start)

	fmt.Printf("Num of times we can generate %v samples from GBM and calc. opt price via FFT for each: %.1f per second.\n",
		numMCSamples, float64(timingN)/elapsed.Seconds())
}
