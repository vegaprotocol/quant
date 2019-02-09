package misc

import (
	"errors"
	"math"
	"math/cmplx"
)

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
