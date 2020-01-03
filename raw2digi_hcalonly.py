import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
import os

process = cms.Process('HCALgpu',eras.Run2_2018)


#-----------------------------------
# Standard CMSSW Imports/Definitions
#-----------------------------------

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = '101X_dataRun2_Prompt_v11'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
# this is with PU (50)
        'file:/data/patatrack/dalfonso/data/2018/Run2018A_HLTPhysics_316457/A4B5F486-5A59-E811-80A4-02163E01A052.root'
### this is the vary high PU
#        'file:/data/patatrack/dalfonso/data/2018/Run2018E_HLTPhysics_325308/FB454F42-97B6-DC4B-88FF-0063C79B9F6C.root'
### this is a standard scenario
#        'file:/data/patatrack/dalfonso/data/2018/Run2018B_HLTPhysics_319300/D6C0583D-5881-E811-9EB8-FA163EAFECF2.root'
    )
)


# load data using the DAQ source
import sys, os, inspect
sys.path.append(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
process.load('sourceFromRaw_cff')


#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange('293765:264-293765:9999')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )



#-----------
# Log output
#-----------
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cout.threshold = "DEBUG"
process.MessageLogger.cerr.threshold = 'DEBUG'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.debugModules = cms.untracked.vstring("*")

process.options = cms.untracked.PSet(
#    numberOfThreads = cms.untracked.uint32(8),
#    numberOfStreams = cms.untracked.uint32(8),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)


process.Out = cms.OutputModule(
        "PoolOutputModule",
        fileName = cms.untracked.string("test.root")
)

#process.Out = cms.OutputModule("AsciiOutputModule",
#                              outputCommands = cms.untracked.vstring(
#       'keep *_hbheprereco_*_*',
#       'keep *_hbheprerecogpu_*_*',
#    ),
#    verbosity = cms.untracked.uint32(0)
#)

#-----------------------------------------
# CMSSW/Hcal non-DQM Related Module import
#-----------------------------------------

process.load("RecoLocalCalo.Configuration.hcalLocalReco_cff")
process.load("EventFilter.HcalRawToDigi.HcalRawToDigi_cfi")
process.load("RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi")

process.hcalDigis.InputLabel = cms.InputTag('rawDataCollector')

#-----------------------------------------
# Temporary customization (things not implemented on the GPU)
#-----------------------------------------

# temporary commenting off -- TimeSlew (M0,MAHI) and PhaseContainment(M0) as not implemented in GPU
process.hbheprereco.algorithm.__setattr__('applyTimeSlew', cms.bool(False))

## the one below is taken directly from the DB, need to comment the line on the CPU code, regard M0
process.hbheprereco.algorithm.__setattr__('correctForPhaseContainment', cms.bool(False))

## do always 8 pulse
process.hbheprereco.algorithm.__setattr__('chiSqSwitch', cms.double(-1))

#-----------------------------------------
# Final Custmization for Run3
#-----------------------------------------

# we will not run arrivla Time at HLT
process.hbheprereco.algorithm.__setattr__('calculateArrivalTime',cms.bool(False))

## the one below is not perfect because of the multiplication for the 1.0f 
process.hbheprereco.algorithm.__setattr__('applyLegacyHBMCorrection', cms.bool(False))

# we only runMahi at HLT
process.hbheprereco.algorithm.__setattr__('useM3',cms.bool(False))

# we will not have the HPD noise flags in Run3, as will be all siPM
process.hbheprereco.setLegacyFlagsQIE8 = cms.bool(False)
process.hbheprereco.setNegativeFlagsQIE8 = cms.bool(False)
process.hbheprereco.setNoiseFlagsQIE8 = cms.bool(False)
process.hbheprereco.setPulseShapeFlagsQIE8 = cms.bool(False)

#-----------------------------------------
#
#-----------------------------------------

# add M0 only for comparison
#process.hbheprerecoM0 = process.hbheprereco.clone()
#process.hbheprerecoM0.algorithm.__setattr__('useMahi',cms.bool(False))

# add analyzer
#process.load("ComparisonPlots.HCALGPUAnalyzer.ConfFile_fig")
#process.comparisonPlots = cms.EDAnalyzer('HCALGPUAnalyzer')
process.TFileService = cms.Service('TFileService', fileName = cms.string('test_v0.root') )


process.hltGetRaw = cms.EDAnalyzer( "HLTGetRaw",
    RawDataCollection = cms.InputTag( "rawDataCollector" )
)

process.input = cms.Path( process.hltGetRaw )

process.finalize = cms.EndPath(process.Out)

process.digiPath = cms.Path(
    process.hcalDigis
)

process.recoPath = cms.Path(
    process.hbheprereco
#    *process.hbheprerecoM0
    *process.hbheprerecogpu
)

process.schedule = cms.Schedule(
#    process.input,
    process.digiPath,
    process.recoPath,
    process.finalize
    )

#dumpFile  = open("dump.py", "w")
#dumpFile.write(process.dumpPython())
#dumpFile.close()
