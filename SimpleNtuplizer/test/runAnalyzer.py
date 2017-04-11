import FWCore.ParameterSet.Config as cms

from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process("EenTest")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32( 1000 )

process.load('Configuration/EventContent/EventContent_cff')
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("RecoEcal.EgammaCoreTools.EcalNextToDeadChannelESProducer_cff")

process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc'   , '')
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 

readFiles.extend([
        '/store/mc/RunIISpring16DR80/DoubleElectron_FlatPt-300To6500/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/00000/0A240112-4406-E611-9A29-20CF3027A5CE.root',
        #'/store/mc/RunIISpring16DR80/DoubleElectron_FlatPt-1To300/AODSIM/PUFlat0to50_80X_mcRun2_asymptotic_2016_v3-v1/00000/0054673A-0306-E611-A949-0025905C2CD0.root'
        #'/store/user/rcoelhol/GluGluHToGG_M-125_13TeV_powheg_pythia8/crab_HggClusters/170327_193908/0000/skimEGMobjects_fromRAW_1.root'
        #'/store/user/rcoelhol/GluGluHToGG_M-125_13TeV_powheg_pythia8/crab_HggClusters_v2/170331_114027/0000/skimEGMobjects_fromRAW_101.root',
])

process.source = cms.Source(
    "PoolSource",
    fileNames = readFiles,
    secondaryFileNames = secFiles
    )

########################################
# Define the analyzer
########################################

process.TFileService = cms.Service("TFileService",
      fileName = cms.string("tree.root"),
      closeFileFast = cms.untracked.bool(True)
  )

process.een_analyzer = cms.EDAnalyzer(
    'SimpleNtuplizer',
    vertices            = cms.InputTag("offlinePrimaryVertices", ""),
    electrons           = cms.InputTag("gedGsfElectrons", ""),
    photons             = cms.InputTag("gedPhotons", ""),
    superClustersEB     = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
    superClustersEE     = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALEndcapWithPreshower"),
    pfLabel             = cms.InputTag("particleFlowClusterECALMatchedToPhotons"),
    rho                 = cms.InputTag("fixedGridRhoFastjetAll", ""),
    genparticles        = cms.InputTag("genParticles", ""),
    PUInfoInputTag      = cms.InputTag("addPileupInfo", ""),
    genEvtInfoInputTag  = cms.InputTag("generator", ""),
    ecalrechitsEB       = cms.InputTag("reducedEcalRecHitsEB"),
    ecalrechitsEE       = cms.InputTag("reducedEcalRecHitsEE"),
    ebSrFlagCollection  = cms.InputTag("ecalDigis"),
    eeSrFlagCollection  = cms.InputTag("ecalDigis"),
    
    doElectronTree      = cms.bool(True),
    doPhotonTree        = cms.bool(True),
    doSuperClusterTree  = cms.bool(True),
    doPFClusterTree     = cms.bool(False),
    doVertex            = cms.bool(True),
    doTagAndProbe       = cms.bool(True),
    saveUnmatched       = cms.bool(False)


    )
    
process.p = cms.Path(process.een_analyzer)
