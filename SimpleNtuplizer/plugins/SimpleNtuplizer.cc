// -*- C++ -*-
//
// Package:    SimpleNtuplizer
// Class:      SimpleNtuplizer
// 
/**\class SimpleNtuplizer SimpleNtuplizer.cc SimpleNtuplizer/plugins/SimpleNtuplizer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


#include "SimpleNtuplizer.h"


////to set the status of gen particle in PF cluster code in a 16 bit integer
void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}

//######################################
//# Constructors and Destructor
//######################################

// Constructor
SimpleNtuplizer::SimpleNtuplizer(const edm::ParameterSet& iConfig):
  // All tokens given in the python config!
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfoInputTag"))),
  genEvtInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfoInputTag"))),
  electronToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  photonToken_(consumes<reco::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
  superClustersEBToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClustersEB"))),
  superClustersEEToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClustersEE"))),
  ecalRecHitEBToken_(consumes<edm::SortedCollection<EcalRecHit>>(iConfig.getParameter<edm::InputTag>("ecalrechitsEB"))),
  ecalRecHitEEToken_(consumes<edm::SortedCollection<EcalRecHit>>(iConfig.getParameter<edm::InputTag>("ecalrechitsEE"))),
  pfLabel_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfLabel"))),
  pspfLabel_(consumes<reco::PFCluster::EEtoPSAssociation>(iConfig.getParameter<edm::InputTag>("pfLabel"))),
  genpfLabel_(consumes<edm::ValueMap<reco::GenParticleRef> >(iConfig.getParameter<edm::InputTag>("pfLabel"))),
  clusSizeLabel_(consumes<edm::ValueMap<int> >(iConfig.getParameter<edm::InputTag>("pfLabel"))),
  ps1Label_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ps1Label"))),
  ps2Label_(consumes<edm::ValueMap<float> >(iConfig.getParameter<edm::InputTag>("ps2Label"))),
  // SRP collections
  ebSrFlagToken_(consumes<EBSrFlagCollection>(iConfig.getParameter<edm::InputTag>("ebSrFlagCollection"))),
  eeSrFlagToken_(consumes<EESrFlagCollection>(iConfig.getParameter<edm::InputTag>("eeSrFlagCollection")))
{
  doElectronTree = iConfig.getParameter<bool>("doElectronTree");
  doPhotonTree = iConfig.getParameter<bool>("doPhotonTree");
  doSuperClusterTree = iConfig.getParameter<bool>("doSuperClusterTree");
  saveUnmatched = iConfig.getParameter<bool>("saveUnmatched");

  doPFTree = iConfig.getParameter<bool>("doPFClusterTree");
  doVertex = iConfig.getParameter<bool>("doVertex");
  doTagAndProbe = iConfig.getParameter<bool>("doTagAndProbe");
  isData = iConfig.getParameter<bool>("isData");

  std::cout << ">>>> Inside SimpleNtuplizer::constructor" << std::endl;

  edm::Service<TFileService> fs;

  eventTree_ = fs->make<TTree> ("EventTree", "Per event data");

  eventTree_->Branch("NtupID_", &NtupID_);

  eventTree_->Branch("eventNumber", &eventNumber_);
  eventTree_->Branch("luminosityBlock", &luminosityBlock_);
  eventTree_->Branch("run", &run_);
  eventTree_->Branch("weight", &weight_);
  eventTree_->Branch("trueNumInteractions", &trueNumInteractions_);	
  if(doVertex) eventTree_->Branch("nPV", &nPV_);
  eventTree_->Branch("nElectrons", &nElectrons_);
  eventTree_->Branch("nElectronsMatched", &nElectronsMatched_);
  eventTree_->Branch("nPhotons", &nPhotons_);
  eventTree_->Branch("nPhotonsMatched", &nPhotonsMatched_);
  eventTree_->Branch("nClusters", &nClusters_);
  eventTree_->Branch("nClustersMatched", &nClustersMatched_);





  if (doElectronTree) {
    electronTree_ = fs->make<TTree> ("ElectronTree", "Electron data");

    electronTree_->Branch("eventNumber", &eventNumber_e);
    electronTree_->Branch("luminosityBlock", &luminosityBlock_e);
    electronTree_->Branch("run", &run_e);
    electronTree_->Branch("isMatched", &isMatched_e);

    electronTree_->Branch("eta", &eta_e);
    electronTree_->Branch("phi", &phi_e);
    electronTree_->Branch("pt", &pt_e);

    electronTree_->Branch("rawEnergy", &rawEnergy_e);
    electronTree_->Branch("preshowerEnergy", &preshowerEnergy_e);
    electronTree_->Branch("preshowerEnergyPlane1", &preshowerEnergyPlane1_e);
    electronTree_->Branch("preshowerEnergyPlane2", &preshowerEnergyPlane2_e);

    electronTree_->Branch("scEta", &scEta_e);
    electronTree_->Branch("scPhi", &scPhi_e);
    electronTree_->Branch("etaWidth", &etaWidth_e);
    electronTree_->Branch("phiWidth", &phiWidth_e);
    electronTree_->Branch("delEtaSeed", &delEtaSeed_e);
    electronTree_->Branch("delPhiSeed", &delPhiSeed_e);

    electronTree_->Branch("isEB", &isEB_e);
    electronTree_->Branch("isEE", &isEE_e);
    electronTree_->Branch("isGap", &isGap_e);
    electronTree_->Branch("isEBEEGap", &isEBEEGap_e);
    electronTree_->Branch("isEBGap", &isEBGap_e);
    electronTree_->Branch("isEBEtaGap", &isEBEtaGap_e);
    electronTree_->Branch("isEBPhiGap", &isEBPhiGap_e);
    electronTree_->Branch("isEEGap", &isEEGap_e);
    electronTree_->Branch("isEEDeeGap", &isEEDeeGap_e);
    electronTree_->Branch("isEERingGap", &isEERingGap_e);

    electronTree_->Branch("rhoValue", &rhoValue_e);

    electronTree_->Branch("seedEnergy", &seedEnergy_e);
    electronTree_->Branch("numberOfClusters", &numberOfClusters_e);

    electronTree_->Branch("clusterRawEnergy", &clusterRawEnergy_e);
    electronTree_->Branch("clusterPhi", &clusterPhi_e);
    electronTree_->Branch("clusterEta", &clusterEta_e);
    electronTree_->Branch("clusterDPhiToSeed", &clusterDPhiToSeed_e);
    electronTree_->Branch("clusterDEtaToSeed", &clusterDEtaToSeed_e);
    electronTree_->Branch("iEtaSeed", &iEtaSeed_e);
    electronTree_->Branch("iPhiSeed", &iPhiSeed_e);
    electronTree_->Branch("iEtaCoordinate", &iEtaCoordinate_e);
    electronTree_->Branch("iPhiCoordinate", &iPhiCoordinate_e);
    electronTree_->Branch("cryEtaCoordinate", &cryEtaCoordinate_e);
    electronTree_->Branch("cryPhiCoordinate", &cryPhiCoordinate_e);
    electronTree_->Branch("iXSeed", &iXSeed_e);
    electronTree_->Branch("iYSeed", &iYSeed_e);
    electronTree_->Branch("iXCoordinate", &iXCoordinate_e);
    electronTree_->Branch("iYCoordinate", &iYCoordinate_e);
    electronTree_->Branch("cryXCoordinate", &cryXCoordinate_e);
    electronTree_->Branch("cryYCoordinate", &cryYCoordinate_e);
    electronTree_->Branch("iEtaMod5", &iEtaMod5_e);
    electronTree_->Branch("iPhiMod2", &iPhiMod2_e);
    electronTree_->Branch("iEtaMod20", &iEtaMod20_e);
    electronTree_->Branch("iPhiMod20", &iPhiMod20_e);
    electronTree_->Branch("r9", &r9_e);
    electronTree_->Branch("sigmaIetaIeta", &sigmaIetaIeta_e);
    electronTree_->Branch("sigmaIetaIphi", &sigmaIetaIphi_e);
    electronTree_->Branch("sigmaIphiIphi", &sigmaIphiIphi_e);
    electronTree_->Branch("e5x5", &e5x5_e);
    electronTree_->Branch("e3x3", &e3x3_e);
    electronTree_->Branch("eMax", &eMax_e);
    electronTree_->Branch("e2nd", &e2nd_e);
    electronTree_->Branch("eTop", &eTop_e);
    electronTree_->Branch("eBottom", &eBottom_e);
    electronTree_->Branch("eLeft", &eLeft_e);
    electronTree_->Branch("eRight", &eRight_e);
    electronTree_->Branch("e2x5Max", &e2x5Max_e);
    electronTree_->Branch("e2x5Left", &e2x5Left_e);
    electronTree_->Branch("e2x5Right", &e2x5Right_e);
    electronTree_->Branch("e2x5Top", &e2x5Top_e);
    electronTree_->Branch("e2x5Bottom", &e2x5Bottom_e);
    electronTree_->Branch("full5x5_r9", &full5x5_r9_e);
    electronTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_e);
    electronTree_->Branch("full5x5_sigmaIetaIphi", &full5x5_sigmaIetaIphi_e);
    electronTree_->Branch("full5x5_sigmaIphiIphi", &full5x5_sigmaIphiIphi_e);
    electronTree_->Branch("full5x5_e5x5", &full5x5_e5x5_e);
    electronTree_->Branch("full5x5_e3x3", &full5x5_e3x3_e);
    electronTree_->Branch("full5x5_eMax", &full5x5_eMax_e);
    electronTree_->Branch("full5x5_e2nd", &full5x5_e2nd_e);
    electronTree_->Branch("full5x5_eTop", &full5x5_eTop_e);
    electronTree_->Branch("full5x5_eBottom", &full5x5_eBottom_e);
    electronTree_->Branch("full5x5_eLeft", &full5x5_eLeft_e);
    electronTree_->Branch("full5x5_eRight", &full5x5_eRight_e);
    electronTree_->Branch("full5x5_e2x5Max", &full5x5_e2x5Max_e);
    electronTree_->Branch("full5x5_e2x5Left", &full5x5_e2x5Left_e);
    electronTree_->Branch("full5x5_e2x5Right", &full5x5_e2x5Right_e);
    electronTree_->Branch("full5x5_e2x5Top", &full5x5_e2x5Top_e);
    electronTree_->Branch("full5x5_e2x5Bottom", &full5x5_e2x5Bottom_e);

    electronTree_->Branch("hadronicOverEm", &hadronicOverEm_e);
    electronTree_->Branch("hadronic1OverEm", &hadronic1OverEm_e);
    electronTree_->Branch("hadronic2OverEm", &hadronic2OverEm_e);
    electronTree_->Branch("hadronicOverEmBC", &hadronicOverEmBC_e);
    electronTree_->Branch("hadronic1OverEmBC", &hadronic1OverEmBC_e);
    electronTree_->Branch("hadronic2OverEmBC", &hadronic2OverEmBC_e);
    electronTree_->Branch("nrSaturatedCrysIn5x5", &nrSaturatedCrysIn5x5_e);
    electronTree_->Branch("nextToDead", &nextToDead_e);

    electronTree_->Branch("corrEnergy74X", &corrEnergy74X_e);
    electronTree_->Branch("corrEnergy74XError", &corrEnergy74XError_e);


    electronTree_->Branch("trkMomentum", &trkMomentum_e);
    electronTree_->Branch("trkMomentumError", &trkMomentumError_e);
    electronTree_->Branch("trkMomentumErrorEGM", &trkMomentumErrorEGM_e);
    electronTree_->Branch("trkMomentumErrorCorrected", &trkMomentumErrorCorrected_e);
    electronTree_->Branch("trkMomentumRelError", &trkMomentumRelError_e);
    electronTree_->Branch("trkEta", &trkEta_e);
    electronTree_->Branch("trkPhi", &trkPhi_e);
    electronTree_->Branch("eOverPuncorr", &eOverPuncorr_e);
    electronTree_->Branch("ecalDriven", &ecalDriven_e);
    electronTree_->Branch("ecalDrivenSeed", &ecalDrivenSeed_e);
    electronTree_->Branch("trackerDrivenSeed", &trackerDrivenSeed_e);
    electronTree_->Branch("classification", &classification_e);
    electronTree_->Branch("fbrem", &fbrem_e);
    electronTree_->Branch("gsfchi2", &gsfchi2_e);
    electronTree_->Branch("gsfndof", &gsfndof_e);
    electronTree_->Branch("gsfnhits", &gsfnhits_e);
	
    electronTree_->Branch("genMatchdR", &genMatchdR_e);
    electronTree_->Branch("genMatchdE", &genMatchdE_e);
    electronTree_->Branch("genMatchdRdE", &genMatchdRdE_e);
    electronTree_->Branch("genPt", &genPt_e);
    electronTree_->Branch("genPhi", &genPhi_e);
    electronTree_->Branch("genEta", &genEta_e);
    electronTree_->Branch("genMass", &genMass_e);
    electronTree_->Branch("genEnergy", &genEnergy_e);
    electronTree_->Branch("genBornEnergy", &genBornEnergy_e);
    electronTree_->Branch("genPdgId", &genPdgId_e);
    electronTree_->Branch("genStatus", &genStatus_e);
    electronTree_->Branch("tp_mll", &tp_mll);
    electronTree_->Branch("tp_ptll", &tp_ptll);
    electronTree_->Branch("tp_tagpt", &tp_tagpt);
    electronTree_->Branch("tp_tageta", &tp_tageta);
    electronTree_->Branch("tp_tagphi", &tp_tagphi);
  }

  if (doSuperClusterTree) {
    superClusterTree_ = fs->make<TTree> ("ClusterTree", "Cluster data");

    superClusterTree_->Branch("eventNumber", &eventNumber_c);
    superClusterTree_->Branch("luminosityBlock", &luminosityBlock_c);
    superClusterTree_->Branch("run", &run_c);
    superClusterTree_->Branch("isMatched", &isMatched_c);

    superClusterTree_->Branch("eta", &eta_c);
    superClusterTree_->Branch("phi", &phi_c);
    superClusterTree_->Branch("pt", &pt_c);

    superClusterTree_->Branch("rawEnergy", &rawEnergy_c);
    superClusterTree_->Branch("preshowerEnergy", &preshowerEnergy_c);
    superClusterTree_->Branch("preshowerEnergyPlane1", &preshowerEnergyPlane1_c);
    superClusterTree_->Branch("preshowerEnergyPlane2", &preshowerEnergyPlane2_c);

    superClusterTree_->Branch("scEta", &scEta_c);
    superClusterTree_->Branch("scPhi", &scPhi_c);
    superClusterTree_->Branch("etaWidth", &etaWidth_c);
    superClusterTree_->Branch("phiWidth", &phiWidth_c);
    superClusterTree_->Branch("delEtaSeed", &delEtaSeed_c);
    superClusterTree_->Branch("delPhiSeed", &delPhiSeed_c);

    superClusterTree_->Branch("isEB", &isEB_c);
    superClusterTree_->Branch("isEE", &isEE_c);
    superClusterTree_->Branch("isGap", &isGap_c);
    superClusterTree_->Branch("isEBEEGap", &isEBEEGap_c);
    superClusterTree_->Branch("isEBGap", &isEBGap_c);
    superClusterTree_->Branch("isEBEtaGap", &isEBEtaGap_c);
    superClusterTree_->Branch("isEBPhiGap", &isEBPhiGap_c);
    superClusterTree_->Branch("isEEGap", &isEEGap_c);
    superClusterTree_->Branch("isEEDeeGap", &isEEDeeGap_c);
    superClusterTree_->Branch("isEERingGap", &isEERingGap_c);

    superClusterTree_->Branch("rhoValue", &rhoValue_c);

    superClusterTree_->Branch("seedEnergy", &seedEnergy_c);
    superClusterTree_->Branch("numberOfClusters", &numberOfClusters_c);

    superClusterTree_->Branch("clusterRawEnergy", &clusterRawEnergy_c);
    superClusterTree_->Branch("clusterPhi", &clusterPhi_c);
    superClusterTree_->Branch("clusterEta", &clusterEta_c);
    superClusterTree_->Branch("clusterDPhiToSeed", &clusterDPhiToSeed_c);
    superClusterTree_->Branch("clusterDEtaToSeed", &clusterDEtaToSeed_c);
    superClusterTree_->Branch("iEtaSeed", &iEtaSeed_c);
    superClusterTree_->Branch("iPhiSeed", &iPhiSeed_c);
    superClusterTree_->Branch("iEtaCoordinate", &iEtaCoordinate_c);
    superClusterTree_->Branch("iPhiCoordinate", &iPhiCoordinate_c);
    superClusterTree_->Branch("cryEtaCoordinate", &cryEtaCoordinate_c);
    superClusterTree_->Branch("cryPhiCoordinate", &cryPhiCoordinate_c);
    superClusterTree_->Branch("iXSeed", &iXSeed_c);
    superClusterTree_->Branch("iYSeed", &iYSeed_c);
    superClusterTree_->Branch("iXCoordinate", &iXCoordinate_c);
    superClusterTree_->Branch("iYCoordinate", &iYCoordinate_c);
    superClusterTree_->Branch("cryXCoordinate", &cryXCoordinate_c);
    superClusterTree_->Branch("cryYCoordinate", &cryYCoordinate_c);
    superClusterTree_->Branch("iEtaMod5", &iEtaMod5_c);
    superClusterTree_->Branch("iPhiMod2", &iPhiMod2_c);
    superClusterTree_->Branch("iEtaMod20", &iEtaMod20_c);
    superClusterTree_->Branch("iPhiMod20", &iPhiMod20_c);
    superClusterTree_->Branch("r9", &r9_c);
    superClusterTree_->Branch("sigmaIetaIeta", &sigmaIetaIeta_c);
    superClusterTree_->Branch("sigmaIetaIphi", &sigmaIetaIphi_c);
    superClusterTree_->Branch("sigmaIphiIphi", &sigmaIphiIphi_c);
    superClusterTree_->Branch("e5x5", &e5x5_c);
    superClusterTree_->Branch("e3x3", &e3x3_c);
    superClusterTree_->Branch("eMax", &eMax_c);
    superClusterTree_->Branch("e2nd", &e2nd_c);
    superClusterTree_->Branch("eTop", &eTop_c);
    superClusterTree_->Branch("eBottom", &eBottom_c);
    superClusterTree_->Branch("eLeft", &eLeft_c);
    superClusterTree_->Branch("eRight", &eRight_c);
    superClusterTree_->Branch("e2x5Max", &e2x5Max_c);
    superClusterTree_->Branch("e2x5Left", &e2x5Left_c);
    superClusterTree_->Branch("e2x5Right", &e2x5Right_c);
    superClusterTree_->Branch("e2x5Top", &e2x5Top_c);
    superClusterTree_->Branch("e2x5Bottom", &e2x5Bottom_c);
    superClusterTree_->Branch("full5x5_r9", &full5x5_r9_c);
    superClusterTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_c);
    superClusterTree_->Branch("full5x5_sigmaIetaIphi", &full5x5_sigmaIetaIphi_c);
    superClusterTree_->Branch("full5x5_sigmaIphiIphi", &full5x5_sigmaIphiIphi_c);
    superClusterTree_->Branch("full5x5_e5x5", &full5x5_e5x5_c);
    superClusterTree_->Branch("full5x5_e3x3", &full5x5_e3x3_c);
    superClusterTree_->Branch("full5x5_eMax", &full5x5_eMax_c);
    superClusterTree_->Branch("full5x5_e2nd", &full5x5_e2nd_c);
    superClusterTree_->Branch("full5x5_eTop", &full5x5_eTop_c);
    superClusterTree_->Branch("full5x5_eBottom", &full5x5_eBottom_c);
    superClusterTree_->Branch("full5x5_eLeft", &full5x5_eLeft_c);
    superClusterTree_->Branch("full5x5_eRight", &full5x5_eRight_c);
    superClusterTree_->Branch("full5x5_e2x5Max", &full5x5_e2x5Max_c);
    superClusterTree_->Branch("full5x5_e2x5Left", &full5x5_e2x5Left_c);
    superClusterTree_->Branch("full5x5_e2x5Right", &full5x5_e2x5Right_c);
    superClusterTree_->Branch("full5x5_e2x5Top", &full5x5_e2x5Top_c);
    superClusterTree_->Branch("full5x5_e2x5Bottom", &full5x5_e2x5Bottom_c);

    superClusterTree_->Branch("nrSaturatedCrysIn5x5", &nrSaturatedCrysIn5x5_c);
    superClusterTree_->Branch("nextToDead", &nextToDead_c);

    superClusterTree_->Branch("corrEnergy74X", &corrEnergy74X_c);
    superClusterTree_->Branch("corrEnergy74XError", &corrEnergy74XError_c);

    superClusterTree_->Branch("genMatchdR", &genMatchdR_c);
    superClusterTree_->Branch("genMatchdE", &genMatchdE_c);
    superClusterTree_->Branch("genMatchdRdE", &genMatchdRdE_c);
    superClusterTree_->Branch("genPt", &genPt_c);
    superClusterTree_->Branch("genPhi", &genPhi_c);
    superClusterTree_->Branch("genEta", &genEta_c);
    superClusterTree_->Branch("genMass", &genMass_c);
    superClusterTree_->Branch("genEnergy", &genEnergy_c);
    superClusterTree_->Branch("genBornEnergy", &genBornEnergy_c);
    superClusterTree_->Branch("genPdgId", &genPdgId_c);
    superClusterTree_->Branch("genStatus", &genStatus_c);
    superClusterTree_->Branch("tp_mll", &tp_mll);
    superClusterTree_->Branch("tp_ptll", &tp_ptll);
    superClusterTree_->Branch("tp_tagpt", &tp_tagpt);
    superClusterTree_->Branch("tp_tageta", &tp_tageta);
    superClusterTree_->Branch("tp_tagphi", &tp_tagphi);
  }

  if (doPhotonTree) {
    photonTree_ = fs->make<TTree> ("PhotonTree", "Photon data");

    photonTree_->Branch("eventNumber", &eventNumber_p);
    photonTree_->Branch("luminosityBlock", &luminosityBlock_p);
    photonTree_->Branch("run", &run_p);
    photonTree_->Branch("isMatched", &isMatched_p);

    photonTree_->Branch("eta", &eta_p);
    photonTree_->Branch("phi", &phi_p);
    photonTree_->Branch("pt", &pt_p);

    photonTree_->Branch("rawEnergy", &rawEnergy_p);
    photonTree_->Branch("preshowerEnergy", &preshowerEnergy_p);
    photonTree_->Branch("preshowerEnergyPlane1", &preshowerEnergyPlane1_p);
    photonTree_->Branch("preshowerEnergyPlane2", &preshowerEnergyPlane2_p);

    photonTree_->Branch("scEta", &scEta_p);
    photonTree_->Branch("scPhi", &scPhi_p);
    photonTree_->Branch("etaWidth", &etaWidth_p);
    photonTree_->Branch("phiWidth", &phiWidth_p);
    photonTree_->Branch("delEtaSeed", &delEtaSeed_p);
    photonTree_->Branch("delPhiSeed", &delPhiSeed_p);

    photonTree_->Branch("isEB", &isEB_p);
    photonTree_->Branch("isEE", &isEE_p);
    photonTree_->Branch("isGap", &isGap_p);
    photonTree_->Branch("isEBEEGap", &isEBEEGap_p);
    photonTree_->Branch("isEBGap", &isEBGap_p);
    photonTree_->Branch("isEBEtaGap", &isEBEtaGap_p);
    photonTree_->Branch("isEBPhiGap", &isEBPhiGap_p);
    photonTree_->Branch("isEEGap", &isEEGap_p);
    photonTree_->Branch("isEEDeeGap", &isEEDeeGap_p);
    photonTree_->Branch("isEERingGap", &isEERingGap_p);

    photonTree_->Branch("rhoValue", &rhoValue_p);

    photonTree_->Branch("seedEnergy", &seedEnergy_p);
    photonTree_->Branch("numberOfClusters", &numberOfClusters_p);

    photonTree_->Branch("clusterRawEnergy", &clusterRawEnergy_p);
    photonTree_->Branch("clusterPhi", &clusterPhi_p);
    photonTree_->Branch("clusterEta", &clusterEta_p);
    photonTree_->Branch("clusterDPhiToSeed", &clusterDPhiToSeed_p);
    photonTree_->Branch("clusterDEtaToSeed", &clusterDEtaToSeed_p);
    photonTree_->Branch("iEtaSeed", &iEtaSeed_p);
    photonTree_->Branch("iPhiSeed", &iPhiSeed_p);
    photonTree_->Branch("iEtaCoordinate", &iEtaCoordinate_p);
    photonTree_->Branch("iPhiCoordinate", &iPhiCoordinate_p);
    photonTree_->Branch("cryEtaCoordinate", &cryEtaCoordinate_p);
    photonTree_->Branch("cryPhiCoordinate", &cryPhiCoordinate_p);
    photonTree_->Branch("iXSeed", &iXSeed_p);
    photonTree_->Branch("iYSeed", &iYSeed_p);
    photonTree_->Branch("iXCoordinate", &iXCoordinate_p);
    photonTree_->Branch("iYCoordinate", &iYCoordinate_p);
    photonTree_->Branch("cryXCoordinate", &cryXCoordinate_p);
    photonTree_->Branch("cryYCoordinate", &cryYCoordinate_p);
    photonTree_->Branch("iEtaMod5", &iEtaMod5_p);
    photonTree_->Branch("iPhiMod2", &iPhiMod2_p);
    photonTree_->Branch("iEtaMod20", &iEtaMod20_p);
    photonTree_->Branch("iPhiMod20", &iPhiMod20_p);
    photonTree_->Branch("r9", &r9_p);
    photonTree_->Branch("sigmaIetaIeta", &sigmaIetaIeta_p);
    photonTree_->Branch("sigmaIetaIphi", &sigmaIetaIphi_p);
    photonTree_->Branch("sigmaIphiIphi", &sigmaIphiIphi_p);
    photonTree_->Branch("e5x5", &e5x5_p);
    photonTree_->Branch("e3x3", &e3x3_p);
    photonTree_->Branch("eMax", &eMax_p);
    photonTree_->Branch("e2nd", &e2nd_p);
    photonTree_->Branch("eTop", &eTop_p);
    photonTree_->Branch("eBottom", &eBottom_p);
    photonTree_->Branch("eLeft", &eLeft_p);
    photonTree_->Branch("eRight", &eRight_p);
    photonTree_->Branch("e2x5Max", &e2x5Max_p);
    photonTree_->Branch("e2x5Left", &e2x5Left_p);
    photonTree_->Branch("e2x5Right", &e2x5Right_p);
    photonTree_->Branch("e2x5Top", &e2x5Top_p);
    photonTree_->Branch("e2x5Bottom", &e2x5Bottom_p);
    photonTree_->Branch("full5x5_r9", &full5x5_r9_p);
    photonTree_->Branch("full5x5_sigmaIetaIeta", &full5x5_sigmaIetaIeta_p);
    photonTree_->Branch("full5x5_sigmaIetaIphi", &full5x5_sigmaIetaIphi_p);
    photonTree_->Branch("full5x5_sigmaIphiIphi", &full5x5_sigmaIphiIphi_p);
    photonTree_->Branch("full5x5_e5x5", &full5x5_e5x5_p);
    photonTree_->Branch("full5x5_e3x3", &full5x5_e3x3_p);
    photonTree_->Branch("full5x5_eMax", &full5x5_eMax_p);
    photonTree_->Branch("full5x5_e2nd", &full5x5_e2nd_p);
    photonTree_->Branch("full5x5_eTop", &full5x5_eTop_p);
    photonTree_->Branch("full5x5_eBottom", &full5x5_eBottom_p);
    photonTree_->Branch("full5x5_eLeft", &full5x5_eLeft_p);
    photonTree_->Branch("full5x5_eRight", &full5x5_eRight_p);
    photonTree_->Branch("full5x5_e2x5Max", &full5x5_e2x5Max_p);
    photonTree_->Branch("full5x5_e2x5Left", &full5x5_e2x5Left_p);
    photonTree_->Branch("full5x5_e2x5Right", &full5x5_e2x5Right_p);
    photonTree_->Branch("full5x5_e2x5Top", &full5x5_e2x5Top_p);
    photonTree_->Branch("full5x5_e2x5Bottom", &full5x5_e2x5Bottom_p);

    photonTree_->Branch("hadronicOverEm", &hadronicOverEm_p);
    photonTree_->Branch("hadronic1OverEm", &hadronic1OverEm_p);
    photonTree_->Branch("hadronic2OverEm", &hadronic2OverEm_p);
    photonTree_->Branch("hadronicOverEmBC", &hadronicOverEmBC_p);
    photonTree_->Branch("hadronic1OverEmBC", &hadronic1OverEmBC_p);
    photonTree_->Branch("hadronic2OverEmBC", &hadronic2OverEmBC_p);
    photonTree_->Branch("nrSaturatedCrysIn5x5", &nrSaturatedCrysIn5x5_p);
    photonTree_->Branch("nextToDead", &nextToDead_p);

    photonTree_->Branch("corrEnergy74X", &corrEnergy74X_p);
    photonTree_->Branch("corrEnergy74XError", &corrEnergy74XError_p);

    photonTree_->Branch("isConverted", &isConverted_p); 

    photonTree_->Branch("genMatchdR", &genMatchdR_p);
    photonTree_->Branch("genMatchdE", &genMatchdE_p);
    photonTree_->Branch("genMatchdRdE", &genMatchdRdE_p);
    photonTree_->Branch("genPt", &genPt_p);
    photonTree_->Branch("genPhi", &genPhi_p);
    photonTree_->Branch("genEta", &genEta_p);
    photonTree_->Branch("genMass", &genMass_p);
    photonTree_->Branch("genEnergy", &genEnergy_p);
    photonTree_->Branch("genBornEnergy", &genBornEnergy_p);
    photonTree_->Branch("genPdgId", &genPdgId_p);
    photonTree_->Branch("genStatus", &genStatus_p);
    photonTree_->Branch("tp_mll", &tp_mll);
    photonTree_->Branch("tp_ptll", &tp_ptll);
    photonTree_->Branch("tp_tagpt", &tp_tagpt);
    photonTree_->Branch("tp_tageta", &tp_tageta);
    photonTree_->Branch("tp_tagphi", &tp_tagphi);
  }


  ////PFtree
  if(doPFTree){

    pfTree_    = fs->make<TTree>("PfTree", "PF Cluster tree");

    pfTree_->Branch("nClus",           &nClus_pf);    
    pfTree_->Branch("clusrawE",        &clusrawE_pf);
    pfTree_->Branch("cluscorrE",       &cluscorrE_pf);
    pfTree_->Branch("clusPt",          &clusPt_pf);
    pfTree_->Branch("clusEta",         &clusEta_pf);
    pfTree_->Branch("clusRho",         &clusRho_pf);
    pfTree_->Branch("clusPhi",         &clusPhi_pf);
    pfTree_->Branch("clusLayer",       &clusLayer_pf);
    pfTree_->Branch("clusSize",        &clusSize_pf);
    pfTree_->Branch("clusIetaIx",      &clusIetaIx_pf);
    pfTree_->Branch("clusIphiIy",      &clusIphiIy_pf);
    pfTree_->Branch("clusPS1",         &clusPS1_pf);
    pfTree_->Branch("clusPS2",         &clusPS2_pf);
    pfTree_->Branch("clusFlag",        &clusFlag_pf);
    pfTree_->Branch("rho",             &rho_pf);
    pfTree_->Branch("nvtx",            &nvtx_pf);
    pfTree_->Branch("genEnergy",       &genEnergy_pf);
    pfTree_->Branch("genPt",           &genPt_pf);
    pfTree_->Branch("genEta",          &genEta_pf);
    pfTree_->Branch("genPhi",          &genPhi_pf);
    pfTree_->Branch("genStatusFlag",          &genStatusFlag_pf);

    pfTree_->Branch("ietamod20",         &ietamod20_pf);
    pfTree_->Branch("iphimod20",         &iphimod20_pf);
    pfTree_->Branch("nhits",         &nhits_pf);
    pfTree_->Branch("tgtvar",         &tgtvar_pf);
    pfTree_->Branch("nlgtgtvar",         &nlgtgtvar_pf);

  }


}



//######################################
//# Member functions
//######################################

// =====================================
// analyze - The method that is executed on every event

void SimpleNtuplizer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ){
    
  using namespace std;
  using namespace edm;
  using namespace reco;

  //######################################
  //# Get all the collections
  //######################################

  // Get vertex collection
  edm::Handle<reco::VertexCollection> vertices;
  if(doVertex) iEvent.getByToken(vtxToken_, vertices);

  // Get electron collection
  edm::Handle<reco::GsfElectronCollection> electrons;  // For AODSIM
  if (doElectronTree)
    iEvent.getByToken(electronToken_, electrons);

  // Get photon collection
  edm::Handle<reco::PhotonCollection> photons;
  if (doPhotonTree)
    iEvent.getByToken(photonToken_, photons);

  // Get clusters  
  edm::Handle<reco::SuperClusterCollection> superClustersEB;    
  edm::Handle<reco::SuperClusterCollection> superClustersEE;    
  if (doSuperClusterTree) {
    iEvent.getByToken( superClustersEBToken_, superClustersEB);
    iEvent.getByToken( superClustersEEToken_, superClustersEE);
  }
  
  if (doElectronTree || doPhotonTree || doSuperClusterTree) {
    iEvent.getByToken( ecalRecHitEBToken_, ecalRecHitsEB_ );
    iEvent.getByToken( ecalRecHitEEToken_, ecalRecHitsEE_ );
  }

  if (!isData) {
    iEvent.getByToken( genParticleToken_, genParticles_ );
    iEvent.getByToken( PUInfoToken_,      puInfoH_ );
    iEvent.getByToken( genEvtInfoToken_,  genEvtInfo_ );      
  }

  edm::ESHandle<CaloGeometry> pGeometry;
  iSetup.get<CaloGeometryRecord>().get(pGeometry);
  geometry_ = pGeometry.product();
    
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  topology_ = pTopology.product();



  //######################################
  //# Event specific quantities (not used in regression)
  //######################################

  // Central event counter (specific to this output tree)
  NtupID_++;

  // Event specific variables
  eventNumber_     = iEvent.id().event();
  luminosityBlock_ = iEvent.id().luminosityBlock();
  run_             = iEvent.id().run();

  if (!isData) {
    weight_ = genEvtInfo_->weight();
    for (std::vector<PileupSummaryInfo>::const_iterator puinfo = puInfoH_->begin(); puinfo != puInfoH_->end(); ++puinfo) {
      if (puinfo->getBunchCrossing() == 0) {
	trueNumInteractions_ = puinfo->getTrueNumInteractions();
	break;
      }
    }
  }
  
  // Determine number of primary vertices
  if(doVertex){
    if (vertices->empty()) nPV_ = 0;
    else nPV_ = vertices->size();
  }

  //######################################
  //# Analyze electrons and photons
  //######################################    

  // Loop over electrons

  if (doElectronTree) {
    nElectrons_ = 0;
    nElectronsMatched_ = 0;
    eventNumber_e     = iEvent.id().event();
    luminosityBlock_e = iEvent.id().luminosityBlock();
    run_e             = iEvent.id().run();
    for (const auto &electron : *electrons) {
      if (doTagAndProbe) 
	if (!findTag(electron, iEvent, iSetup)) return;
      setElectronVariables(electron, iEvent, iSetup);	
    }
  }

  // Loop over photons
  if (doPhotonTree) {
    nPhotons_         = 0;
    nPhotonsMatched_  = 0;
    eventNumber_p     = iEvent.id().event();
    luminosityBlock_p = iEvent.id().luminosityBlock();
    run_p             = iEvent.id().run();

    for (const auto &photon : *photons) {
      if (doTagAndProbe) 
	if (!findTag(photon, iEvent, iSetup)) return;
      setPhotonVariables( photon, iEvent, iSetup );
    }
  }


  if (doSuperClusterTree) {
    nClusters_ = 0;
    nClustersMatched_ = 0;
    eventNumber_c     = iEvent.id().event();
    luminosityBlock_c = iEvent.id().luminosityBlock();
    run_c             = iEvent.id().run();
    for (const auto &superCluster : *superClustersEB) {
      setSuperClusterVariables( superCluster, iEvent, iSetup, true );
    }
    for (const auto &superCluster : *superClustersEE) {
      if (doTagAndProbe) { // Some gymnastics needed here
	math::XYZPoint v(0, 0, 0);
	math::XYZVector p = superCluster.energy() * (superCluster.position() - v).unit();
	reco::RecoEcalCandidate recoSuperCluster(1, reco::Candidate::LorentzVector(p.x(), p.y(), p.z(), sqrt(p.mag2())), v, 11, 1);
	if (!findTag(recoSuperCluster, iEvent, iSetup)) return;
      }
      setSuperClusterVariables( superCluster, iEvent, iSetup, false );
    }
  }
    

  ///doPFTree
  if(doPFTree){
    setPFVariables(iEvent, iSetup );
  }

  
  // Fill in the event specific variables
  eventTree_->Fill();

}

//######################################
//# Necessary functions and settings for framework
//######################################

// ------------ method called once each job just before starting event loop  ------------
void SimpleNtuplizer::beginJob() {
  //std::cout << ">>>> Inside SimpleNtuplizer::beginJob" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void SimpleNtuplizer::endJob() {
  //std::cout << ">>>> Inside SimpleNtuplizer::endJob" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleNtuplizer);

