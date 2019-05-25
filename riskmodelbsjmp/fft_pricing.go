package riskmodelbsjmp

import (
	"math"
	"math/cmplx"

	"code.vegaprotocol.io/quant/misc"
)

const i1 complex128 = complex(0, 1)

// CallPriceFftHalfAxis given characteristic function of underlying at time T this method will calculate the
// call price using the FFT for S = 1 and "all" strikes
// This method is described in Vega Technical paper: Margins and Credit Risk
func CallPriceFftHalfAxis(zeta func(float64) complex128, r, T float64) ([]float64, []float64, error) {
	const log2ofN int = 17
	minK := math.Pow(10, -12)
	N := int(math.Pow(2, float64(log2ofN)))
	minSmallK := math.Log(minK)

	kGrid, f, err := misc.FFTOfEvenInRealPartFn(minSmallK, N, zeta)
	if err != nil {
		return make([]float64, 0), make([]float64, 0), err
	}

	// post process
	const log2ofLittleN int = 9 // how many we will drop at the out of money end
	littleN := int(math.Pow(2, float64(log2ofLittleN)))

	strikesGrid := make([]float64, len(kGrid)-littleN)
	pricesGrid := make([]float64, len(kGrid)-littleN)
	for i := 0; i < len(strikesGrid); i++ {
		j := i + littleN
		strikesGrid[i] = math.Exp(kGrid[j])
		pricesGrid[i] = real(f[j]) + math.Max(1.0-math.Exp(kGrid[j]-r*T), 0.0)
	}
	return strikesGrid, pricesGrid, err
}

// CallPriceFftFullAxis given characteristic function of underlying at time T this method will calculate the
// call price using the FFT for S = 1 and "all" strikes
// This method is described e.g. in
// Cont, R. and Tankov, P. - Financial Modelling with Jump Processes, Section 11.1.3
func CallPriceFftFullAxis(zeta func(float64) complex128, r, T float64) ([]float64, []float64, error) {

	const A float64 = 10000
	const log2ofN int = 17

	N := int(math.Pow(2, float64(log2ofN)))
	const L float64 = 4.0 / 2.0 / math.Pi

	un := make([]float64, N)
	kDeltaL := make([]float64, N)
	fftIn := make([]complex128, N)

	Delta := A / float64(N-1)
	for i := 0; i < N; i++ {
		vk := -A/2 + float64(i)*Delta
		wk := 1.0
		un[i] = 2 * math.Pi * (float64(i)/(float64(N)*Delta) - L)
		kDeltaL[i] = float64(i) * Delta * 2 * math.Pi * L
		if i == 0 {
			wk = 0.5
		}
		if i == N-1 {
			wk = 0.5
		}
		fftIn[i] = complex(wk, 0) * zeta(vk) * cmplx.Exp(1i*complex(kDeltaL[i], 0))
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
