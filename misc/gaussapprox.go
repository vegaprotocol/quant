package misc

import "math"

// GaussDensity returns the density of N(0,1) r.v.
func GaussDensity(x float64) float64 {
	A := 1.0 / math.Sqrt(2.0*math.Pi)
	return A * math.Exp(-x*x*0.5)
}

// ApproxGaussCdf returns a (fast) approximation to the distribution of N(0,1)
func ApproxGaussCdf(x float64) float64 {
	const a1 float64 = 0.4361836
	const a2 float64 = -0.1201676
	const a3 float64 = 0.9372980

	var k = 1.0 / (1.0 + (0.33267 * x))
	if x >= 0.0 {
		return 1.0 - GaussDensity(x)*(a1*k+(a2*k*k)+(a3*k*k*k))
	}
	return 1.0 - ApproxGaussCdf(-x)
}
