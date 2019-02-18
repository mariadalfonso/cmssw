import FWCore.ParameterSet.Config as cms

trackingMaterialAnalyser = cms.EDAnalyzer("TrackingMaterialAnalyser",
    MaterialAccounting      = cms.InputTag("trackingMaterialProducer"),
    SplitMode               = cms.string("NearestLayer"),
    SkipBeforeFirstDetector = cms.bool(False),
    SkipAfterLastDetector   = cms.bool(True),
    SaveSummaryPlot         = cms.bool(True),
    SaveDetailedPlots       = cms.bool(False),
    SaveParameters          = cms.bool(True),
    SaveXML                 = cms.bool(True),
# to derive the following list:
# cat ../data/trackingMaterialGroups_ForPhaseII.xml | grep TrackingMaterialGroup | sed -e 's/\s*//' | cut -d ' ' -f 3 | tr '=' ' ' | cut -d ' ' -f 2 | tr -d '"' | sed -e 's/\(.*\)/"\1",/'
    Groups = cms.vstring(
        "HGCalEESensitive",
        "HGCalHESiliconSensitive",
        "HGCalHEScintillatorSensitive"
    )
)
