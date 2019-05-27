package fftpricing

import (
	"math"
	"math/cmplx"

	"code.vegaprotocol.io/quant/bsformula"
	"code.vegaprotocol.io/quant/misc"
)

const i1 complex128 = complex(0, 1)

// PricingModelParamsBSJmp collect the parameters of Black-Scholes model.
// Here R is the risk-free interest rate, Sigma is volatiliy
// Gamma is jump intensity (jump inter-arrival times are exponential with this param)
// finally the jump sizes are normally distributed N(a,b^2),
// JmpMeanA = a, variance = JmpStddevB^2 = b^2
type PricingModelParamsBSJmp struct {
	R          float64 // interest rate
	Sigma      float64 // volatility of diffusion part
	Gamma      float64 // jump intensity (jump inter-arrival times are exponential with this param)
	JmpMeanA   float64 // jump mean
	JmpStddevB float64 // jump stddev
}

// CreatePricingParamsForBSJmp creates an instance of PricingModelParamsBSJmp
func CreatePricingParamsForBSJmp(r, sigma, gamma, jmpMeanA, jmpStddevB float64) PricingModelParamsBSJmp {
	var paramsJmp PricingModelParamsBSJmp
	paramsJmp.R = r
	paramsJmp.Sigma = sigma
	paramsJmp.Gamma = gamma
	paramsJmp.JmpMeanA = jmpMeanA
	paramsJmp.JmpStddevB = jmpStddevB
	return paramsJmp
}

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

// CallPriceFftFullAxis given characteristic function zeta for the distribution
// of the underlying at time T this method will calculate the
// call price using the FFT for S = 1 and "all" strikes
// it needs the correct "integration corrector function"
// This method is described e.g. in
// Cont, R. and Tankov, P. - Financial Modelling with Jump Processes, Section 11.1.3
func CallPriceFftFullAxis(zeta func(float64) complex128, correction func(float64) float64) ([]float64, []float64, error) {

	const A float64 = 10000
	const log2ofN int = 15

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
	const log2ofLittleN int = 12 // how many we will drop at the out of money end
	littleN := int(math.Pow(2, float64(log2ofLittleN)))
	strikes := make([]float64, N-littleN)
	prices := make([]float64, N-littleN)
	for i := 0; i < len(prices); i++ {
		j := i + littleN
		strikes[i] = math.Exp(un[j])
		//prices[i] = real(fk[j]) + math.Max(1.0-math.Exp(un[j]-r*T), 0)
		prices[i] = real(fk[j]) + correction(un[j])
	}
	return strikes, prices, err
}

// r interest rate
// sigma volatility of diffusion part
// gamma jump intensity (jump inter-arrival times are exponential with this param)
// jmpMeanA   float64 jump mean
// jmpStddevB float64 jump stddev of normally dist. jumps
func jumpDiffusionCharacteristicFn(u complex128, T, r, sigma, gamma, jmpMeanA, jmpStddevB float64) complex128 {
	alpha := math.Exp(jmpMeanA+0.5*jmpStddevB*jmpStddevB) - 1.0
	firstExponential := cmplx.Exp(complex(0, (r-0.5*sigma*sigma-alpha*gamma)*T)*u - complex(0.5*sigma*sigma*T, 0)*u*u)
	secondExponential := complex(1, 0)
	if gamma > 0 {
		secondExponential = cmplx.Exp(complex(gamma*T, 0) * (misc.CharacteristicOfGaussian(u, r, sigma) - 1))
	}
	return firstExponential * secondExponential
}

// JumpDiffCallPricesFftFullAxis calculates risk-neutral call price for an option
// with current asset price S = 1 for "all" strikes
func JumpDiffCallPricesFftFullAxis(T float64, p PricingModelParamsBSJmp) ([]float64, []float64, error) {
	r := p.R
	sigma := p.Sigma
	//halfSigmaSquareT := complex(0.5*sigma*sigma*T, 0)

	// characteristic of jump diff part
	phiOne := func(u complex128) complex128 {
		return jumpDiffusionCharacteristicFn(u, T, r, sigma, p.Gamma, p.JmpMeanA, p.JmpStddevB)
	}
	// characteristic of normal r.v. with variance given by sigma^2 * T
	phiTwo := func(u complex128) complex128 {
		return jumpDiffusionCharacteristicFn(u, T, r, sigma, 0, 0, 0)
	}

	// this is the relevant characteristic function
	zeta := func(v float64) complex128 {
		vC := complex(v, 0)
		return cmplx.Exp(complex(0, v*r*T)) * (phiOne(vC-i1) - phiTwo(vC-i1)) / (i1 * vC * (1 + i1*vC))
	}

	// correction function
	correction := func(u float64) float64 {
		// the input u here is the log-moneyness
		K := math.Exp(u)
		return bsformula.BSCallPrice(1.0, K, r, sigma, T)
	}
	return CallPriceFftFullAxis(zeta, correction)
}
