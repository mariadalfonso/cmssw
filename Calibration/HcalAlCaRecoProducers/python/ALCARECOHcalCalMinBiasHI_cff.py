import FWCore.ParameterSet.Config as cms

#------------------------------------------------
#AlCaReco filtering for HCAL minbias:
#------------------------------------------------

from Calibration.HcalAlCaRecoProducers.ALCARECOHcalCalMinBiasNoise_cff import *
from Calibration.HcalAlCaRecoProducers.ALCARECOHcalCalMinBias_cff import *

import HLTrigger.HLTfilters.hltHighLevel_cfi
hcalminbiasHLT =  HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
#    HLTPaths = ['HLT_HcalPhiSym'],
    eventSetupPathsKey='HcalCalMinBiasHI',
    throw = False #dont throw except on unknown path name 
)

import RecoLocalCalo.HcalRecProducers.HBHEPhase1Reconstructor_cfi


seqALCARECOHcalCalMinBiasDigi = cms.Sequence(hcalminbiasHLT*hcalDigiAlCaMB*gtDigisAlCaMB)
seqALCARECOHcalCalMinBias = cms.Sequence(hbherecoNoise*hfrecoNoise*hfrecoMBNZS*horecoNoise)


_phase1_seqALCARECOHcalCalMinBias = seqALCARECOHcalCalMinBias.copy()
_phase1_seqALCARECOHcalCalMinBias.insert(0,hfprerecoMBNZS)
_phase1_seqALCARECOHcalCalMinBias.insert(0,hfprerecoNoise)

from Configuration.Eras.Modifier_run2_HF_2017_cff import run2_HF_2017
run2_HF_2017.toReplaceWith( seqALCARECOHcalCalMinBias, _phase1_seqALCARECOHcalCalMinBias )
run2_HF_2017.toReplaceWith( hfrecoNoise, phase1_hfrecoNoise )
run2_HF_2017.toReplaceWith( hfrecoMBNZS, phase1_hfrecoMBNZS )

_plan1_seqALCARECOHcalCalMinBias = _phase1_seqALCARECOHcalCalMinBias.copy()
hbheprerecoNoise = hbherecoNoise.clone()
hbheprerecoMBNZS = hbherecoMBNZS.clone()
_plan1_seqALCARECOHcalCalMinBias.insert(0,hbheprerecoNoise)
_plan1_seqALCARECOHcalCalMinBias.insert(0,hbheprerecoMBNZS)
from Configuration.Eras.Modifier_run2_HEPlan1_2017_cff import run2_HEPlan1_2017
run2_HEPlan1_2017.toReplaceWith(hbherecoNoise, hbheplan1Noise)
run2_HEPlan1_2017.toReplaceWith(hbherecoMBNZS, hbheplan1MBNZS)
run2_HEPlan1_2017.toReplaceWith(seqALCARECOHcalCalMinBias, _plan1_seqALCARECOHcalCalMinBias)
