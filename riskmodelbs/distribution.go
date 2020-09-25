package riskmodelbs

import (
	"math"

	"code.vegaprotocol.io/quant/interfaces"

	"gonum.org/v1/gonum/stat/distuv"
)

const probabilityTolerance = 1e-3

// GetProbabilityDistribution returns the log normal distribution correcposning to the suplied model parameters, the current stock price S and time horizon tau.
func (modelParams ModelParamsBS) GetProbabilityDistribution(S, tau float64) interfaces.AnalyticalDistribution {
	// See: http://www-2.rotman.utoronto.ca/~hull/TechnicalNotes/TechnicalNote2.pdf
	m := math.Log(S) + (modelParams.Mu-0.5*modelParams.Sigma*modelParams.Sigma)*tau
	stdDev := modelParams.Sigma * math.Sqrt(tau)

	return &distuv.LogNormal{Mu: m, Sigma: stdDev}
}

// GetProbabilityTolerance specifies the probability tolerance alphaModel that the model supports. It shouldn't be used for any calculations involving alpha < alphaModel or alpha > 1-alphaModel
func (modelParams ModelParamsBS) GetProbabilityTolerance() (alphaModel float64) {
	alphaModel = probabilityTolerance
	return alphaModel
}
