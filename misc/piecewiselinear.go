package misc

import "math"

// Function is a piecewise-linear 1-dimensional function
type Function struct {
	X []float64
	Y []float64
}

// Area returns the definite integral of the function on its domain X.
//
// Time complexity: O(N), where N is the number of points.
// Space complexity: O(1)
func (f Function) Area() (area float64) {
	X, Y := f.X, f.Y
	for i := 1; i < len(X); i++ {
		area += (X[i] - X[i-1]) * (Y[i] + Y[i-1]) / 2
	}
	return area
}

// AreaUpTo returns the definite integral of the function on its domain X intersected with [-Inf, x].
//
// Time complexity: O(N), where N is the number of points.
// Space complexity: O(1)
func (f Function) AreaUpTo(x float64) (area float64) {
	X, Y := f.X, f.Y
	for i := 1; i < len(X); i++ {
		dX := X[i] - X[i-1]
		if x < X[i] {
			if x >= X[i-1] {
				dxX := x - X[i-1]
				w := dxX / dX
				y := (1-w)*Y[i-1] + w*Y[i]
				area += dxX * (y + Y[i-1]) / 2
			}
			return area
		}
		area += dX * (Y[i] + Y[i-1]) / 2
	}
	return area
}

// At returns the function's value at the given point.
// Outside its domain X, the function is constant at 0.
//
// The function's X and Y slices are expected to be the same legnth. The length property is _not_ verified.
// The function's X slice is expected to be sorted in ascending order. The sortedness property is _not_ verified.
//
// Time complexity: O(log(N)), where N is the number of points.
// Space complexity: O(1)
func (f Function) At(x float64) float64 {
	X, Y := f.X, f.Y
	i, j := 0, len(X)
	for i < j {
		h := int(uint(i+j) >> 1)
		if X[h] < x {
			i = h + 1
		} else {
			j = h
		}
	}
	switch i {
	case 0, len(X):
		return 0
	}
	w := (x - X[i-1]) / (X[i] - X[i-1])
	return (1-w)*Y[i-1] + w*Y[i]
}

// Span generates `nPoints` equidistant points spanning [min,max]
func Span(min, max float64, nPoints int) []float64 {
	X := make([]float64, nPoints)
	min, max = math.Min(max, min), math.Max(max, min)
	d := max - min
	for i := range X {
		X[i] = min + d*(float64(i)/float64(nPoints-1))
	}
	return X
}
