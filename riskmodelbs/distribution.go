package riskmodelbs

import (
	"math"

	"code.vegaprotocol.io/quant/interfaces"

	"gonum.org/v1/gonum/stat/distuv"
)

const probabilityTolerance = 1e-3

func (modelParams ModelParamsBS) GetProbabilityDistribution(S, tau float64) interfaces.AnalyticalDistribution {
	mean := S * math.Exp(modelParams.Mu*tau)
	stdDev := S * math.Exp(modelParams.Mu*tau) * math.Sqrt(math.Exp(modelParams.Sigma*modelParams.Sigma*tau)-1)
	return &distuv.LogNormal{Mu: mean, Sigma: stdDev}
}

func (modelParams ModelParamsBS) GetProbabilityTolerance() float64 {
	return probabilityTolerance
}
