import FWCore.ParameterSet.Config as cms

# helper fuctions
from HLTrigger.Configuration.common import *

# add one customisation function per PR
# - put the PR number into the name of the function
# - add a short comment
# for example:

# CCCTF tuning
# def customiseFor12718(process):
#     for pset in process._Process__psets.values():
#         if hasattr(pset,'ComponentType'):
#             if (pset.ComponentType == 'CkfBaseTrajectoryFilter'):
#                 if not hasattr(pset,'minGoodStripCharge'):
#                     pset.minGoodStripCharge = cms.PSet(refToPSet_ = cms.string('HLTSiStripClusterChargeCutNone'))
#     return process

def customiseFor21821(process):
    from RecoLocalCalo.HcalRecProducers.HBHEPhase1Reconstructor_cfi import hbheprereco
    for producer in producers_by_type(process, "HBHEPhase1Reconstructor"):
        producer.algorithm.ts4Max = cms.vdouble(100., 20000., 30000)
        del producer.algorithm.pedestalUpperLimit
        del producer.algorithm.pedSigmaHPD
        del producer.algorithm.pedSigmaSiPM
        del producer.algorithm.noiseHPD
        del producer.algorithm.noiseSiPM

    for producer in producers_by_type(process, "HcalHitReconstructor"):
        if hasattr(producer,"puCorrMethod"):
            del producer.puCorrMethod

    return process

def customiseFor21664_forMahiOn(process):
    for producer in producers_by_type(process, "HBHEPhase1Reconstructor"):
        producer.algorithm.useMahi   = cms.bool(True)
        producer.algorithm.useM2     = cms.bool(False)
        producer.algorithm.useM3     = cms.bool(False)
    return process

def customiseFor21664_forMahiOnM2only(process):
    for producer in producers_by_type(process, "HBHEPhase1Reconstructor"):
      if (producer.algorithm.useM2 == cms.bool(True)):
        producer.algorithm.useMahi   = cms.bool(True)
        producer.algorithm.useM2     = cms.bool(False)
        producer.algorithm.useM3     = cms.bool(False)
    return process

# CMSSW version specific customizations
def customizeHLTforCMSSW(process, menuType="GRun"):

    # add call to action function in proper order: newest last!
    # process = customiseFor12718(process)

    process = customiseFor21821(process)

    return process
