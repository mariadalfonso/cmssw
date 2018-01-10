import FWCore.ParameterSet.Config as cms
from RecoParticleFlow.PFClusterProducer.particleFlowCaloResolution_cfi import _timeResolutionHCAL

particleFlowRecHitHBHE = cms.EDProducer("PFRecHitProducer",
    navigator = cms.PSet(
            name = cms.string("PFRecHitHCALNavigator"),
            sigmaCut = cms.double(4.0),
            timeResolutionCalc = _timeResolutionHCAL
    ),
    producers = cms.VPSet(
           cms.PSet(
             name = cms.string("PFHBHERecHitCreator"),
             src  = cms.InputTag("hbhereco",""),
             qualityTests = cms.VPSet(
                  cms.PSet(
                  name = cms.string("PFRecHitQTestHCALThresholdVsDepth"),
                  cuts = cms.VPSet(
                        cms.PSet(
                            depth=cms.double(1.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalBarrel"),
                            ),
                        cms.PSet(
                            depth=cms.double(2.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalBarrel"),
                            ),
                        cms.PSet(
                            depth=cms.double(3.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalBarrel"),
                            ),
                        cms.PSet(
                            depth=cms.double(1.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalEndcap"),
                            ),
                        cms.PSet(
                            depth=cms.double(2.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalEndcap"),
                            )
                        cms.PSet(
                            depth=cms.double(3.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalEndcap"),
                            )
                        cms.PSet(
                            depth=cms.double(4.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalEndcap"),
                            )
                        cms.PSet(
                            depth=cms.double(5.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalEndcap"),
                            )
                        cms.PSet(
                            depth=cms.double(6.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalEndcap"),
                            )
                        cms.PSet(
                            depth=cms.double(7.),
                            threshold = cms.double(0.8),
                            detector = cms.string("HcalEndcap"),
                            )
                        )
                  ),
                  cms.PSet(
                      name = cms.string("PFRecHitQTestHCALChannel"),
                      maxSeverities      = cms.vint32(11),
                      cleaningThresholds = cms.vdouble(0.0),
                      flags              = cms.vstring('Standard')
                  )
                  

             )
           ),
           
    )

)

