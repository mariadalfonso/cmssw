import FWCore.ParameterSet.Config as cms

#------------------------------------------------
#AlCaReco filtering for HCAL pedestal:
#------------------------------------------------

import HLTrigger.HLTfilters.hltHighLevel_cfi
hcalCalibPedestalHLT =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    eventSetupPathsKey='HcalCalPedestal',
    throw = False #dont throw except on unknown path name 
)

import EventFilter.HcalRawToDigi.HcalCalibTypeFilter_cfi
hcalCalibPedestal = EventFilter.HcalRawToDigi.HcalCalibTypeFilter_cfi.hcalCalibTypeFilter.clone(
    #  InputLabel = cms.string('rawDataCollector'),
    InputLabel = cms.string('hltHcalCalibrationRaw::HLT'),
    #  InputLabel = cms.InputTag("hltEcalCalibrationRaw","","HLT"),
    CalibTypes    = cms.vint32( 1 ),
    FilterSummary = cms.untracked.bool( False )
    )

#add GT digi:
import EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi
gtDigisAlCaPedestal = EventFilter.L1GlobalTriggerRawToDigi.l1GtUnpack_cfi.l1GtUnpack.clone()

import EventFilter.HcalRawToDigi.HcalRawToDigi_cfi
hcalDigiAlCaPedestal = EventFilter.HcalRawToDigi.HcalRawToDigi_cfi.hcalDigis.clone()
hcalDigiAlCaPedestal.InputLabel = cms.InputTag('hltHcalCalibrationRaw')

qie10Digis = EventFilter.HcalRawToDigi.HcalRawToDigi_cfi.hcalDigis.clone()
qie10Digis.InputLabel = cms.InputTag('hltHcalCalibrationRaw')
qie10Digis.FEDs = cms.untracked.vint32(1132)


from Calibration.HcalAlCaRecoProducers.ALCARECOHcalCalMinBiasNoise_cff import *

hbherecoPedestal = hbherecoNoise.clone()
hbherecoPedestal.digiLabelQIE8  = cms.InputTag("hcalDigiAlCaPedestal")
hbherecoPedestal.digiLabelQIE11 = cms.InputTag("hcalDigiAlCaPedestal")

horecoPedestal = horecoNoise.clone()
horecoPedestal.digiLabel = cms.InputTag('hcalDigiAlCaPedestal')

hfrecoPedestal = hfrecoNoise.clone()
hfrecoPedestal.digiLabel = cms.InputTag('hcalDigiAlCaPedestal')

hfprerecoPedestal = hfprerecoNoise.clone()
hfprerecoPedestal.digiLabel = cms.InputTag("hcalDigiAlCaPedestal")

phase1_hfrecoPedestal = phase1_hfrecoNoise.clone()
phase1_hfrecoPedestal.inputLabel = cms.InputTag("hfprerecoPedestal")



seqALCARECOHcalCalPedestal = cms.Sequence(hbherecoPedestal*hfrecoPedestal*horecoPedestal)

seqALCARECOHcalCalPedestalDigi = cms.Sequence(hcalCalibPedestalHLT*
                                              hcalCalibPedestal*
                                              hcalDigiAlCaPedestal*
                                              qie10Digis*
                                              gtDigisAlCaPedestal)

_phase1_seqALCARECOHcalCalPedestal = seqALCARECOHcalCalPedestal.copy()
_phase1_seqALCARECOHcalCalPedestal.insert(0,hfprerecoPedestal)

from Configuration.Eras.Modifier_run2_HF_2017_cff import run2_HF_2017
run2_HF_2017.toReplaceWith( seqALCARECOHcalCalPedestal, _phase1_seqALCARECOHcalCalPedestal )
run2_HF_2017.toReplaceWith( hfrecoPedestal, phase1_hfrecoPedestal )

import RecoLocalCalo.HcalRecProducers.hbheplan1_cfi
hbheplan1Pedestal = RecoLocalCalo.HcalRecProducers.hbheplan1_cfi.hbheplan1.clone(
    hbheInput = cms.InputTag("hbheprerecoPedestal")
)

_plan1_seqALCARECOHcalCalPedestal = _phase1_seqALCARECOHcalCalPedestal.copy()
hbheprerecoPedestal = hbherecoPedestal.clone()
_plan1_seqALCARECOHcalCalPedestal.insert(0,hbheprerecoPedestal)
from Configuration.Eras.Modifier_run2_HEPlan1_2017_cff import run2_HEPlan1_2017
run2_HEPlan1_2017.toReplaceWith(hbherecoPedestal, hbheplan1Pedestal)
run2_HEPlan1_2017.toReplaceWith(seqALCARECOHcalCalPedestal, _plan1_seqALCARECOHcalCalPedestal)
