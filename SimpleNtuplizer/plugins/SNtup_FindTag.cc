#include "SimpleNtuplizer.h"

#include "FWCore/Utilities/interface/RegexMatch.h"
#include "DataFormats/L1TGlobal/interface/GlobalObjectMap.h"
#include "DataFormats/L1TGlobal/interface/GlobalObjectMapRecord.h"
#include "DataFormats/L1TGlobal/interface/GlobalObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include <boost/foreach.hpp>
#include <iostream>

bool SimpleNtuplizer::findTag(
        const reco::RecoCandidate& object,
        const edm::Event& iEvent,
        const edm::EventSetup& iSetup ){

  using namespace std;

  tp_mll = -1;
  tp_ptll = -1;
  tp_tagpt = -1;
  tp_tageta = -1;
  tp_tagphi = -1;
  
  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> electrons;  // For AODSIM
  iEvent.getByToken(electronToken_, electrons);

  // Loop over electrons
  bool foundIt = false;
  for( auto&& electron : *electrons) {

    bool passKinematics = electron.pt() > 40. && fabs(electron.eta()) < 2.5;
    if (!passKinematics) continue;
    
    bool passID = electron.full5x5_sigmaIetaIeta() < 0.01 &&
      electron.deltaEtaSeedClusterTrackAtVtx() < 0.003 &&
      electron.deltaPhiSuperClusterTrackAtVtx() < 0.08 &&
      electron.hcalOverEcal() < 0.01 &&
      (1./electron.correctedEcalEnergy() - 1./electron.trackMomentumAtVtx().rho()) < 0.01 &&
      //electron.gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS) < 1 &&      
      electron.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS) < 1 &&      
      (electron.pfIsolationVariables().sumChargedHadronPt + electron.pfIsolationVariables().sumNeutralHadronEt + electron.pfIsolationVariables().sumPhotonEt) < 0.05*electron.pt();
    if (!passID) continue;

    bool passB2B = cos(reco::deltaPhi(electron.phi(), object.phi())) < -0.99;
    if (!passB2B) continue;
    
    foundIt = true;

    float mass = (electron.p4()+object.p4()).mass();
    float pt = (electron.p4()+object.p4()).pt();
    tp_mll = mass;
    tp_ptll = pt;
    tp_tagpt = electron.pt();
    tp_tageta = electron.eta();
    tp_tagphi = electron.phi();

  }
  return foundIt;
    
}
