package riskmodelbs

import (
	"math"

	"gonum.org/v1/gonum/stat/distuv"
)

const ProbabilityTolerance = 1e-3

func GetDistribution(modelParams ModelParamsBS, S, tau float64) *distuv.LogNormal {
	mean := S * math.Exp(modelParams.Mu*tau)
	stDev := S * S * math.Exp(2*modelParams.Mu*tau) * (math.Exp(modelParams.Sigma*modelParams.Sigma*tau) - 1)
	return &distuv.LogNormal{Mu: mean, Sigma: stDev}
}
