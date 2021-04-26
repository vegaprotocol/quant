package interfaces

// AnalyticalDistribution ...
type AnalyticalDistribution interface {
	Mean() float64
	Variance() float64
	CDF(float64) float64
	Quantile(float64) float64
}

// AnalyticalModel ...
type AnalyticalModel interface {
	GetProbabilityDistribution(S, tau float64) AnalyticalDistribution
	GetProbabilityTolerance() float64
}
