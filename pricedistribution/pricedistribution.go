package pricedistribution

type AnalyticDistribution interface {
	CDF(float64) float64
	Quantile(float64) float64
}

func PriceRange(d AnalyticDistribution, alpha float64) (Smin float64, Smax float64) {
	x := (1 - alpha) / 2
	Smin = d.Quantile(x)
	Smax = d.Quantile(alpha + x)
	return
}

func PriceDistribution(d AnalyticDistribution, bins []float64) (probabilities []float64) {
	probabilities = make([]float64, 0, len(bins)-1)
	pLeft := d.CDF(bins[0])
	for i := 1; i < len(bins); i++ {
		pRight := d.CDF(bins[i])
		pBin := pRight - pLeft
		probabilities = append(probabilities, pBin)
		pLeft = pRight
	}
	return
}
