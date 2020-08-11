package riskmodelbs

import (
	"testing"

	"code.vegaprotocol.io/quant/interfaces"
)

func TestBSIsAnalyticalModel(t *testing.T) {
	const r float64 = 0.0
	const mu float64 = 0.0
	const sigma float64 = 1.5

	const tau = 1.0 / 365 / 24 / 60
	const lambda float64 = 0.01

	var analyticBs interfaces.AnalyticalModel = ModelParamsBS{Mu: 0, R: 0, Sigma: 1.5}

	if analyticBs == nil {
		t.Error("Expeced ModelParamsBS to implement AnalyticalModel interface")
	}
}
