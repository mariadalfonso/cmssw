import FWCore.ParameterSet.Config as cms
import math

L1TTriggerTowerConfig_etaphi = cms.PSet(readMappingFile=cms.bool(False),
                                        minEta=cms.double(1.479),
                                        maxEta=cms.double(3.0),
                                        minPhi=cms.double(-1*math.pi),
                                        maxPhi=cms.double(math.pi),
                                        nBinsEta=cms.int32(18),
                                        nBinsPhi=cms.int32(72),
                                        binsEta=cms.vdouble(),
                                        binsPhi=cms.vdouble())

towerMap2D_parValues = cms.PSet( useLayerWeights = cms.bool(False),
                                 layerWeights = cms.vdouble(),
                                 L1TTriggerTowerConfig = L1TTriggerTowerConfig_etaphi
)

tower_map = cms.PSet( ProcessorName  = cms.string('HGCalTowerMapProcessor'),
                      towermap_parameters = towerMap2D_parValues.clone()
                  )

hgcalTowerMapProducer = cms.EDProducer(
    "HGCalTowerMapProducer",
    InputTriggerSums = cms.InputTag('hgcalConcentratorProducer:HGCalConcentratorProcessorSelection'),
    ProcessorParameters = tower_map.clone()
    )

### update the sequence for the HFNose

L1TTriggerTowerConfigHFNose_etaphi = L1TTriggerTowerConfig_etaphi.clone(
    minEta=cms.double(3.0),
    maxEta=cms.double(4.2)
)

towerMap2DHFNose_parValues = towerMap2D_parValues.clone(
    L1TTriggerTowerConfig = L1TTriggerTowerConfigHFNose_etaphi
)


towerHFNose_map = cms.PSet( ProcessorName  = cms.string('HGCalTowerMapProcessor'),
                      towermap_parameters = towerMap2DHFNose_parValues.clone()
                  )
hgcalTowerMapProducerHFNose = hgcalTowerMapProducer.clone(
    InputTriggerSums = cms.InputTag('hgcalConcentratorProducerHFNose:HGCalConcentratorProcessorSelection'),
    ProcessorParameters = towerHFNose_map.clone()
)
