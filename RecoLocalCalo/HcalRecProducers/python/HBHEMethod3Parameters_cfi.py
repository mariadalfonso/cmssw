import FWCore.ParameterSet.Config as cms

# Configuration parameters for Method 3
m3Parameters = cms.PSet(
    applyTimeSlewM3         = cms.bool(True),
    respCorrM3              = cms.double(1.0)   # This factor is used to align the Method3 with the Method2 response
)
