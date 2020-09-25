package pricedistribution

import (
	"code.vegaprotocol.io/quant/interfaces"
)

// PriceRange returns the minimum and maximum price implied by the supplied distribution and probability level
func PriceRange(d interfaces.AnalyticalDistribution, alpha float64) (minPrice float64, maxPrice float64) {
	x := (1 - alpha) / 2
	minPrice = d.Quantile(x)
	maxPrice = d.Quantile(alpha + x)
	return
}

// PriceDistribution returns the probability implied by the supplied distribution per each interval specified by the bins,
// when less than 2 bins are specifed an empty array gets returned.
func PriceDistribution(d interfaces.AnalyticalDistribution, bins []float64) (probabilities []float64) {
	n := len(bins)
	if n < 2 {
		probabilities = make([]float64, 0, 0)
		return
	}

	probabilities = make([]float64, 0, n-1)
	pLeft := d.CDF(bins[0])
	for _, bin := range bins[1:] {
		pRight := d.CDF(bin)
		pBin := pRight - pLeft
		probabilities = append(probabilities, pBin)
		pLeft = pRight
	}
	return
}

// ProbabilityOfTrading returns a probability of trading implied by the supplied distribution, price and order side (bid/offer). minPrice and maxPrice are used to constrain the distribution
// in case applyMinMax flag is set to true.
func ProbabilityOfTrading(d interfaces.AnalyticalDistribution, price float64, isBid bool, applyMinMax bool, minPrice float64, maxPrice float64) float64 {
	if isBid {
		if applyMinMax && (price <= minPrice || price > maxPrice) {
			return 0
		}
		min := 0.0
		if applyMinMax {
			min = d.CDF(minPrice)
		}
		return d.CDF(price) - min
	}
	if applyMinMax && (price < minPrice || price >= maxPrice) {
		return 0
	}
	max := 1.0
	if applyMinMax {
		max = d.CDF(maxPrice)
	}
	return max - d.CDF(price)

}
