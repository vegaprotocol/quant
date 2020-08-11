package pricedistribution

func PriceRange(invCDF func(float64) float64, alpha float64) (Smin float64, Smax float64) {
	x := (1 - alpha) / 2
	Smin = invCDF(x)
	Smax = invCDF(alpha + x)
	return
}

func PriceDistribution(CDF func(float64) float64, bins []float64) (probabilities []float64) {
	probabilities = make([]float64, 0, len(bins)-1)
	pLeft := CDF(bins[0])
	for i := 1; i < len(bins); i++ {
		pRight := CDF(bins[i])
		pBin := pRight - pLeft
		probabilities = append(probabilities, pBin)
		pLeft = pRight
	}
	return
}
