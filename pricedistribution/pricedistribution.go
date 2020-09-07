package pricedistribution

import (
	"code.vegaprotocol.io/quant/interfaces"
)

func PriceRange(d interfaces.AnalyticalDistribution, alpha float64) (Smin float64, Smax float64) {
	x := (1 - alpha) / 2
	Smin = d.Quantile(x)
	Smax = d.Quantile(alpha + x)
	return
}

func PriceDistribution(d interfaces.AnalyticalDistribution, bins []float64) (probabilities []float64) {
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

func ProbabilityOfTrading(d interfaces.AnalyticalDistribution, price float64, isBid bool, Smin float64, Smax float64) float64 {
	if isBid {
		if price <= Smin || price > Smax {
			return 0
		}
		return d.CDF(price) - d.CDF(Smin)
	}
	if price < Smin || price >= Smax {
		return 0
	}
	cdfSmax := d.CDF(Smax)
	cdfPrice := d.CDF(price)
	return cdfSmax - cdfPrice

}
