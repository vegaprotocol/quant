package misc

import "math"
import "errors"

const h float64 = 1e-6

//FindRootWithoutDerivative returns an approximate solution s to f(x)=0 using Newtons method starting at x0 and such that |f(x)| < maxError
//where f is a function R->R
//Results in error if derivative is too small or number of iterations exceeds maxIter.
func FindRootWithoutDerivative(f func(float64) float64,
	x0 float64, maxIter int, maxError float64) (float64, error) {

	fPrime := func(x float64) float64 {
		return (f(x+h) - f(x-h)) / (2 * h)
	}
	return FindRoot(f, fPrime, x0, maxIter, maxError)
}

//FindRoot returns an approximate solution s to f(x)=0 using Newtons method starting at x0 and such that |f(x)| < maxError
//where f is a function R->R and fPrime is f'
//Results in error if derivative is too small or number of iterations exceeds maxIter.
func FindRoot(f, fPrime func(float64) float64,
	x0 float64, maxIter int, maxError float64) (float64, error) {

	var i int //  = 0
	xn := x0

	for i < maxIter {
		if math.Abs(f(xn)) < maxError {
			return xn, nil
		}
		if math.Abs(fPrime(xn)) < 1e-16 {
			return math.NaN(), errors.New("NewtonsMethod failed - derivative too small")
		}
		xn = xn - f(xn)/fPrime(xn)
		i++
	}

	return math.NaN(), errors.New("NewtonsMethod did not converge")
}
