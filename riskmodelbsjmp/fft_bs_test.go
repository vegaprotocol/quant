package riskmodelbsjmp

import (
	//"fmt"
	"math"
	"math/cmplx"
	"testing"

	"code.vegaprotocol.io/quant/bsformula"
	"code.vegaprotocol.io/quant/misc"
)

// TestFFTMethodForBlackScholes can calculate the price of call option in the Black-Scholes
// model using FFT; so if the FFT code is correct then
// the price should closely match the BS formula
func TestFFTMethodForBlackScholes(t *testing.T) {
	const tol float64 = 1e-5
	const S float64 = 1
	const minK float64 = 0.0000001
	const sigma float64 = 0.5
	const r float64 = 0.0
	const T float64 = 10

	strikesGrid, optVals, err := bsCallPriceFft(r, sigma, T)
	if err != nil {
		t.Errorf(err.Error())
		return
	}

	bsValsFFT := misc.Function{
		X: strikesGrid,
		Y: optVals,
	}
	K := 1.0
	bsAtmPrice := bsformula.BSCallPrice(S, K, r, sigma, T)
	fftAtmPrice := bsValsFFT.At(K)
	if math.Abs(bsAtmPrice-fftAtmPrice) > tol {
		t.Errorf("BS = %g, FFT = %g", bsAtmPrice, fftAtmPrice)
	}
}

const i1 complex128 = complex(0, 1)

func TestFFTMethodForBlackScholes2(t *testing.T) {
	// the very out of money and very in-the-money values are not that accurate
	const tol float64 = 1e-3
	const sigma float64 = 0.5
	const r float64 = 0.0
	const T float64 = 10
	const S float64 = 1

	strikes, optVals, err := bsCallPriceFft(r, sigma, T)
	if err != nil {
		t.Errorf(err.Error())
		return
	}
	N := len(optVals)
	for i := 0; i < N; i++ {
		bsPrice := bsformula.BSCallPrice(S, strikes[i], r, sigma, T)
		fftPrice := optVals[i]
		if math.Abs(bsPrice-fftPrice) > tol {
			t.Errorf("BS = %g, FFT = %g at K = %g", bsPrice, fftPrice, strikes[i])
			return
		}
	}
}

// Fourier transform price for S = 1 and "all" strikes
func bsCallPriceFft(r, sigma, T float64) ([]float64, []float64, error) {

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

	const A float64 = 10000
	const log2ofN int = 15

	N := int(math.Pow(2, float64(log2ofN)))
	const L float64 = 4.0 / 2.0 / math.Pi

	vk := make([]float64, N)
	wk := make([]float64, N) // weights
	un := make([]float64, N)
	kDeltaL := make([]float64, N)
	fftIn := make([]complex128, N)

	Delta := A / float64(N-1)
	for i := 0; i < N; i++ {
		vk[i] = -A/2 + float64(i)*Delta
		wk[i] = 1.0
		un[i] = 2 * math.Pi * (float64(i)/(float64(N)*Delta) - L)
		kDeltaL[i] = float64(i) * Delta * 2 * math.Pi * L
		if i == 0 {
			wk[i] = 0.5
		}
		if i == N-1 {
			wk[i] = 0.5
		}
		fftIn[i] = complex(wk[i], 0) * zeta(vk[i]) * cmplx.Exp(1i*complex(kDeltaL[i], 0))
	}
	fftOut, err := misc.FFT(fftIn)
	if err != nil {
		return make([]float64, 0), make([]float64, 0), err
	}

	fk := make([]complex128, N)
	for i := 0; i < N; i++ {
		fk[i] = complex(A/(float64(N)*2*math.Pi), 0) * cmplx.Exp(1i*complex(un[i]*0.5*A, 0)) * fftOut[i]
	}

	// post process
	const log2ofLittleN int = 9 // how many we will drop at the out of money end
	littleN := int(math.Pow(2, float64(log2ofLittleN)))
	strikes := make([]float64, N-littleN)
	prices := make([]float64, N-littleN)
	for i := 0; i < len(prices); i++ {
		j := i + littleN
		strikes[i] = math.Exp(un[j])
		prices[i] = real(fk[j]) + math.Max(1.0-math.Exp(un[j]-r*T), 0)
	}
	return strikes, prices, err
}
