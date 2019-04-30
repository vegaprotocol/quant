package misc

import (
	"math"
	"math/cmplx"
	"testing"

	"gitlab.com/vega-protocol/quant/bsformula"
)

func TestFFTInputChecks(t *testing.T) {
	yVals := make([]complex128, 0)
	yHatVals, err := fft(yVals)
	if err == nil {
		t.Errorf("FFT should fail for arrays of lenght 0")
	}
	if len(yHatVals) != 0 {
		t.Errorf("returned array must have lenght 0")
	}

	yVals2 := make([]complex128, 11)
	yHatVals2, err2 := fft(yVals2)
	if err2 == nil {
		t.Errorf("FFT should fail for arrays of lenght that are not power of 2")
	}
	if len(yHatVals2) != 0 {
		t.Errorf("returned array must have lenght 0")
	}

}

func TestFourierAndInverse(t *testing.T) {
	const tolerance float64 = 1e-12
	N := int(math.Pow(2.0, 16))
	A := 10.0
	L := 5.0
	xGrid := make([]float64, N)
	yVals := make([]complex128, N)
	Delta := A / float64(N-1)
	for i := 0; i < N; i++ {
		xGrid[i] = 2 * math.Pi * (Delta*float64(i) - L)
		yVals[i] = complex(math.Sin(xGrid[i]), 0)
	}

	yHatVals, errFft := fft(yVals)
	if errFft != nil {
		t.Errorf(errFft.Error())
		return
	}
	yValsAfter, errIfft := ifft(yHatVals)
	if errIfft != nil {
		t.Errorf(errFft.Error())
		return
	}
	for i := 0; i < N; i++ {
		absDiff := cmplx.Abs(yVals[i] - yValsAfter[i])
		if absDiff > tolerance {
			t.Errorf("Diff = %v at index i = %v", absDiff, i)
			return
		}
	}
}

func TestFFTvsPresaved(t *testing.T) {
	const tolerance float64 = 1e-6
	N := 16
	A := 10.0
	L := 5.0
	xGrid := make([]float64, N)
	yVals := make([]complex128, N)
	Delta := A / float64(N-1)
	for i := 0; i < N; i++ {
		xGrid[i] = 2 * math.Pi * (Delta*float64(i) - L)
		yVals[i] = complex(math.Sin(xGrid[i]), 0)
	}

	pythonyHatValsReal := [16]float64{1.26565425e-14, -2.31488188e-02, -1.05066500e-01, -3.02833275e-01, -8.66025404e-01, -5.10345215e+00, 3.56916811e+00, 1.96533263e+00, 1.73205081e+00, 1.96533263e+00, 3.56916811e+00, -5.10345215e+00, -8.66025404e-01, -3.02833275e-01, -1.05066500e-01, -2.31488188e-02}
	pythonyHatValsImag := [16]float64{0.0, 0.11637697, 0.25365297, 0.45322202, 0.8660254, 3.4100177, -1.47839784, -0.39092897, 0., 0.39092897, 1.47839784, -3.4100177, -0.8660254, -0.45322202, -0.25365297, -0.11637697}

	yHatVals, err := fft(yVals)
	if err != nil {
		t.Errorf(err.Error())
	} else {
		for i := 0; i < N; i++ {
			if math.Abs(pythonyHatValsReal[i]-real(yHatVals[i])) > tolerance {
				t.Errorf("Error real part: python FFT %g, misc.FFT %g", pythonyHatValsReal[i], real(yHatVals[i]))
			}
			if math.Abs(pythonyHatValsImag[i]-imag(yHatVals[i])) > tolerance {
				t.Errorf("Error imag part: python FFT %g, misc.FFT %g", pythonyHatValsImag[i], imag(yHatVals[i]))
			}
		}
	}
}

// We can calculate the price of call option in the Black-Scholes
// model using FFT; so if the FFT code is correct then
// the price should closely match the BS formula
func TestFFTMethodForBlackScholes(t *testing.T) {
	const tol float64 = 1e-2
	const S float64 = 100
	const minK float64 = 1
	const sigma float64 = 0.1
	const r float64 = 0.0
	const T float64 = 1

	strikesGrid, optVals, err := BSCallPriceFft(minK, S, r, sigma, T)
	if err != nil {
		t.Errorf(err.Error())
		return
	}

	BSValsFFT := Function{
		X: strikesGrid,
		Y: optVals,
	}
	K := 100.0
	bsAtmPrice := bsformula.BSCallPrice(S, K, r, sigma, T)
	fftAtmPrice := BSValsFFT.At(K)
	if math.Abs(bsAtmPrice-fftAtmPrice) > tol {
		t.Errorf("BS = %g, FFT = %g", bsAtmPrice, fftAtmPrice)
	}
}

const i1 complex128 = complex(0, 1)

func BSCallPriceFft(minK, S, r, sigma, T float64) ([]float64, []float64, error) {
	N := int(math.Pow(2, 14))
	minSmallK := math.Log(minK / S)
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
	gridSmallK, f, err := FourierTransformOfEvenInRealPartFn(minSmallK, N, zeta)
	if err != nil {
		return make([]float64, 0), make([]float64, 0), err
	}
	optVals := make([]float64, N)
	strikesGrid := make([]float64, N)
	for i := 0; i < N; i++ {
		optVals[i] = S * (real(f[i]) + math.Max(1.0-math.Exp(gridSmallK[i]-r*T), 0.0))
		strikesGrid[i] = S * math.Exp(gridSmallK[i])
	}
	return strikesGrid, optVals, err
}
