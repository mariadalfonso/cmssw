// -*- C++ -*-
//
// Package:    SimpleHcalNtuplizer
// Class:      SimpleHcalNtuplizer
// 
/**\class SimpleHcalNtuplizer SimpleHcalNtuplizer.cc SimpleHcalNtuplizer/plugins/SimpleHcalNtuplizer.cc

   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Ilya Kravchenko
//         Created:  Thu, 10 Jul 2014 09:54:13 GMT
//
//


#include "SimpleHcalNtuplizer.h"

/*
////to set the status of gen particle in PF cluster code in a 16 bit integer
void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}
*/

//######################################
//# Constructors and Destructor
//######################################

// Constructor
SimpleHcalNtuplizer::SimpleHcalNtuplizer(const edm::ParameterSet& iConfig):
  // All tokens given in the python config!
  vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genparticles"))),
  //  rhoToken_(consumes<double> (iConfig.getParameter<edm::InputTag>("rho"))),
  PUInfoToken_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("PUInfoInputTag"))),
  genEvtInfoToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfoInputTag"))),
  pfLabel_(consumes<reco::PFClusterCollection>(iConfig.getParameter<edm::InputTag>("pfLabel")))
{

  doSuperClusterTree = false;
  doElectronTree = false;
  doPhotonTree = false;


  std::cout << ">>>> Inside SimpleHcalNtuplizer::constructor" << std::endl;

  edm::Service<TFileService> fs;

  eventTree_ = fs->make<TTree> ("EventTree", "Per event data");

  eventTree_->Branch("NtupID_", &NtupID_);

  eventTree_->Branch("eventNumber", &eventNumber_);
  eventTree_->Branch("luminosityBlock", &luminosityBlock_);
  eventTree_->Branch("run", &run_);
  eventTree_->Branch("weight", &weight_);
  eventTree_->Branch("trueNumInteractions", &trueNumInteractions_);	
  if(doVertex) eventTree_->Branch("nPV", &nPV_);
  //  eventTree_->Branch("nElectrons", &nElectrons_);
  //  eventTree_->Branch("nElectronsMatched", &nElectronsMatched_);
  //  eventTree_->Branch("nPhotons", &nPhotons_);
  //  eventTree_->Branch("nPhotonsMatched", &nPhotonsMatched_);
  eventTree_->Branch("nClusters", &nClusters_);
  eventTree_->Branch("nClustersMatched", &nClustersMatched_);


  ////PFtree
  if(doPFTree){

    pfTree_    = fs->make<TTree>("PfTree", "PF Cluster tree");

    pfTree_->Branch("nClus",           &nClus_pf);    
    pfTree_->Branch("clusrawE",        &clusrawE_pf);
    pfTree_->Branch("cluscorrE",       &cluscorrE_pf);
    pfTree_->Branch("clusPt",          &clusPt_pf);
    pfTree_->Branch("clusEta",         &clusEta_pf);
    pfTree_->Branch("clusHits",        &clusHits_pf);
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
    pfTree_->Branch("genStatusFlag",          &genStatusFlag_pf);

    pfTree_->Branch("genMatchdR", &genMatchdR_e);
    pfTree_->Branch("genMatchdE", &genMatchdE_e);
    pfTree_->Branch("genMatchdRdE", &genMatchdRdE_e);
    pfTree_->Branch("genPt", &genPt_e);
    pfTree_->Branch("genPhi", &genPhi_e);
    pfTree_->Branch("genEta", &genEta_e);
    pfTree_->Branch("genEnergy", &genEnergy_e);
    pfTree_->Branch("genMass", &genMass_e);
    pfTree_->Branch("genBornEnergy", &genBornEnergy_e);
    pfTree_->Branch("genPdgId", &genPdgId_e);
    pfTree_->Branch("genStatus", &genStatus_e);

  }


}



//######################################
//# Member functions
//######################################

// =====================================
// analyze - The method that is executed on every event

void SimpleHcalNtuplizer::analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup ){
    
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

  /*
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
  */    

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
void SimpleHcalNtuplizer::beginJob() {
  //std::cout << ">>>> Inside SimpleHcalNtuplizer::beginJob" << std::endl;
}

// ------------ method called once each job just after ending the event loop  ------------
void SimpleHcalNtuplizer::endJob() {
  //std::cout << ">>>> Inside SimpleHcalNtuplizer::endJob" << std::endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SimpleHcalNtuplizer);

