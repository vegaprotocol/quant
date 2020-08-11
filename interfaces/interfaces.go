package interfaces

type AnalyticalDistribution interface {
	CDF(float64) float64
	Quantile(float64) float64
}

type AnalyticalModel interface {
	GetProbabilityDistribution(S, tau float64) AnalyticalDistribution
	GetProbabilityTolerance() float64
}