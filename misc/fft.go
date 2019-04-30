package misc

import (
	"errors"
	"math"
	"math/cmplx"
)

////////////////////////////////////////////////////////////////////////////////
// Inverse of the Fast Fourier Transform
// for a vector a in an N dimensional space of complex numbers,
// of lenght which is a power of 2,
// this calculates, for all k = 0,1,...,N-1
// X_k = (1/N)\sum_{n=0}^{N-1} a_n exp( 2 \pi i k n / N )
//
// It uses the fft method and the observation that we can get the
// inverse of FFT by conjugating the input and then conjugating the output
////////////////////////////////////////////////////////////////////////////////
func ifft(a []complex128) ([]complex128, error) {
	n := len(a)
	y := make([]complex128, n)

	// to do an inverse transform we first conjugate the input...
	for i := 0; i < n; i++ {
		y[i] = complex(real(a[i]), -1.0*imag(a[i]))
	}
	// ... do the forward fourier transform on this ...
	y, err := fft(y)
	// ... conjugate the output
	for i := 0; i < n; i++ {
		y[i] = complex(real(y[i])/float64(n), -1.0*imag(y[i])/float64(n))
	}
	return y, err
}

////////////////////////////////////////////////////////////////////////////////
// Fast Fourier Transform
// for a vector a in an N dimensional space of complex numbers,
// this calculates, for all k = 0,1,...,N-1
// of lenght which is a power of 2,
// X_k = \sum_{n=0}^{N-1} a_n exp( -2 \pi i k n / N )
//
// However instead of a naive algorithm that would need N^2 operations
// this only needs N log N. The catch is that the input lenght must be power of 2
//
////////////////////////////////////////////////////////////////////////////////
func fft(a []complex128) ([]complex128, error) {
	n := len(a)
	if n <= 0 {
		return make([]complex128, 0), errors.New("Array length must be > 0")
	}
	if (n & (n - 1)) != 0 {
		return make([]complex128, 0), errors.New("Array length must be power of 2")
	}

	x := make([]complex128, n)
	copy(x, a)

	j := 0
	for i := 0; i < n; i++ {
		if i < j {
			x[i], x[j] = x[j], x[i]
		}
		m := n / 2
		for {
			if j < m {
				break
			}
			j = j - m
			m = m / 2
			if m < 2 {
				break
			}
		}
		j = j + m
	}
	kmax := 1
	for {
		if kmax >= n {
			return x, nil
		}
		istep := kmax * 2
		for k := 0; k < kmax; k++ {
			theta := complex(0.0, -1.0*math.Pi*float64(k)/float64(kmax))
			for i := k; i < n; i += istep {
				j := i + kmax
				temp := x[j] * cmplx.Exp(theta)
				x[j] = x[i] - temp
				x[i] = x[i] + temp
			}
		}
		kmax = istep
	}
}

// FourierTransformOfEvenInRealPartFn will do the fourier transform of a given function
// xMin - minimum value in the grid in Fourier space, must be -ve!
// N - number of discretization points, must be power of 2
// f - the function to do fourier transform of, it must be odd in imaginary part
// even in the real part (so we integrate over (0,\infty) only)
func FourierTransformOfEvenInRealPartFn(
	xiMin float64, N int, f func(float64) complex128) ([]float64, []complex128, error) {

	if N <= 0 {
		return make([]float64, 0), make([]complex128, 0), errors.New("N must be > 0")
	}
	if (N & (N - 1)) != 0 {
		return make([]float64, 0), make([]complex128, 0), errors.New("N must be power of 2")
	}
	if xiMin >= 0.0 {
		return make([]float64, 0), make([]complex128, 0), errors.New("xiMin must be strictly negative")
	}

	DeltaXi := (-xiMin - xiMin) / float64(N-1)
	xiGrid := make([]float64, N) // grid in the fourier space
	DeltaX := 2.0 * math.Pi / (float64(N) * DeltaXi)
	xGrid := make([]float64, N)
	a := make([]complex128, N)
	scaleFactor := (1.0 / math.Pi) * DeltaX
	for i := 0; i < N; i++ {
		xiGrid[i] = xiMin + float64(i)*DeltaXi
		xGrid[i] = float64(i+1) * DeltaX
		fTimesExp := f(xGrid[i]) * cmplx.Exp(complex(0.0, xiMin*DeltaX*float64(i)))
		a[i] = complex(scaleFactor*real(fTimesExp), scaleFactor*imag(fTimesExp))
	}
	fHat, err := fft(a)

	return xiGrid, fHat, err
}
