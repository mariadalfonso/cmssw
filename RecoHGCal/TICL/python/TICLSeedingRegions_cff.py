import FWCore.ParameterSet.Config as cms

from RecoHGCal.TICL.ticlSeedingRegionProducer_cfi import ticlSeedingRegionProducer as _ticlSeedingRegionProducer

# SEEDING REGION

ticlSeedingGlobal = _ticlSeedingRegionProducer.clone(
  algoId = 2
)

ticlSeedingTrk = _ticlSeedingRegionProducer.clone(
  algoId = 1
)

ticlSeedingGlobalHFNose = _ticlSeedingRegionProducer.clone(
  algoId = 2
)

ticlSeedingTrkHFNose = _ticlSeedingRegionProducer.clone(
    algoId = 1,
    detector = "HFNose",
    cutTk = cms.string('3. < abs(eta) < 4. && pt > 1. && quality("highPurity") && hitPattern().numberOfLostHits("MISSING_OUTER_HITS") < 5'),
)
