package misc

//"gonum.org/v1/gonum/stat/"
import (
	"gonum.org/v1/gonum/stat"
)

// CalculateControlVariateEstimator does
func CalculateControlVariateEstimator(unknownSamples []float64, knownSamples []float64, knownMean float64) float64 {
	N := len(unknownSamples)
	empMeanOfUnknown := stat.Mean(unknownSamples, nil)
	empMeanOfKnown := stat.Mean(knownSamples, nil)
	bStarNominator := 0.0
	bStarDenominator := 0.0
	for i := 0; i < N; i++ {
		// x refers to the known
		// y refers to the unknown
		xCentred := knownSamples[i] - empMeanOfKnown
		yCentred := unknownSamples[i] - empMeanOfUnknown
		bStarNominator += xCentred*yCentred
		bStarDenominator += xCentred*xCentred
	}
	bStar := bStarNominator/bStarDenominator 

	controlVariateEstimator := 0.0
	for i := 0; i < N; i++ {
		controlVariateEstimator += unknownSamples[i] - bStar*(knownSamples[i] - knownMean)
	}
	controlVariateEstimator = controlVariateEstimator / float64(N)
	return controlVariateEstimator
}
