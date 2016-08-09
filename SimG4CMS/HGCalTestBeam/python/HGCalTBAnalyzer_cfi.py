import FWCore.ParameterSet.Config as cms

HGCalTBAnalyzer = cms.EDAnalyzer("HGCalTBAnalyzer",
                                 DetectorEE   = cms.string('HGCalEESensitive'),
                                 UseEE        = cms.bool(True),
                                 DetectorHE   = cms.string('HGCalHESiliconSensitive'),
                                 UseHE        = cms.bool(False),
                                 DoSimHits    = cms.bool(True),
                                 DoDigis      = cms.bool(True),
                                 DoRecHits    = cms.bool(True),
                                 SampleIndex  = cms.int32(0),
                                 GeneratorSrc = cms.InputTag('generatorSmeared'),
                                 CaloHitSrcEE = cms.string('HGCHitsEE'),
                                 CaloHitSrcHE = cms.string('HGCHitsHEfront'),
                                 DigiSrcEE    = cms.InputTag('mix','HGCDigisEE'),
                                 DigiSrcHE    = cms.InputTag('mix','HGCDigisHEfront'),
                                 RecHitSrcEE  = cms.InputTag('HGCalRecHit','HGCEERecHits'),
                                 RecHitSrcHE  = cms.InputTag('HGCalRecHit','HGCHEFRecHits'),
                                 MakeTree     = cms.untracked.bool(True)
)
