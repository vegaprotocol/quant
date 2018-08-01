package riskmeasures

import (
	"math"

	"gonum.org/v1/gonum/stat/distuv"
)

// LogNormalVaR computes value at risk of LogNormal r.v.
func LogNormalVaR(mu, sigma, alpha float64) float64 {
	return -math.Exp(mu + sigma*distuv.UnitNormal.Quantile(alpha))
}

// NegLogNormalVaR computes value at risk of LogNormal r.v.
func NegativeLogNormalVaR(mu, sigma, alpha float64) float64 {
	return math.Exp(mu + sigma*distuv.UnitNormal.Quantile(1.0-alpha))
}

//LogNormalEs returns the expected shortfall of a lognormal r.v. at given lambda level
func LogNormalEs(mu, sigma, lambd float64) float64 {
	//var x = distuv.UnitNormal.CDF(0)
	var quantileForLambda = distuv.UnitNormal.Quantile(lambd)
	var es = -(1.0 / lambd) * math.Exp(mu+sigma*sigma*0.5) * distuv.UnitNormal.CDF(quantileForLambda-sigma)
	return es
}

//NegativeLogNormal returns the expected shortfall of a lognormal r.v. at given lambda level
func NegativeLogNormalEs(mu, sigma, lambd float64) float64 {
	//var x = distuv.UnitNormal.CDF(0)
	var quantileForOneMinusLambda = distuv.UnitNormal.Quantile(1.0 - lambd)
	var es = (1.0 / lambd) * math.Exp(mu+sigma*sigma*0.5) * (1 - distuv.UnitNormal.CDF(quantileForOneMinusLambda-sigma))
	return es
}
