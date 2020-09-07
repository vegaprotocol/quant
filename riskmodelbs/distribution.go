package riskmodelbs

import (
	"math"

	"code.vegaprotocol.io/quant/interfaces"

	"gonum.org/v1/gonum/stat/distuv"
)

const probabilityTolerance = 1e-3

func (modelParams ModelParamsBS) GetProbabilityDistribution(S, tau float64) interfaces.AnalyticalDistribution {
	// See: http://www-2.rotman.utoronto.ca/~hull/TechnicalNotes/TechnicalNote2.pdf
	m := math.Log(S) + (modelParams.Mu-0.5*modelParams.Sigma*modelParams.Sigma)*tau
	stdDev := modelParams.Sigma * math.Sqrt(tau)

	return &distuv.LogNormal{Mu: m, Sigma: stdDev}
}

func (modelParams ModelParamsBS) GetProbabilityTolerance() float64 {
	return probabilityTolerance
}
