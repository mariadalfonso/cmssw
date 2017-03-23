#include "SimpleNtuplizer.h"

// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setSuperClusterVariables(
					       const reco::SuperCluster& superCluster,
					       const edm::Event& iEvent,
					       const edm::EventSetup& iSetup,
					       bool isEB){

  using namespace std;
  using namespace edm;
  using namespace reco;

  // Increase count of photons in event
  nClusters_++;

  // Try to match to genParticle; quit function if superCluster is not matched
  isMatched_c = matchSuperClusterToGenParticle(superCluster);    
  if (!isMatched_c && !saveUnmatched) return;
  else nClustersMatched_ += isMatched_c ? 1 : 0;

  // Get the right ecalRecHits collection (different for barrel and encap)
  edm::Handle<SortedCollection<EcalRecHit>> ecalRecHits;
  if (isEB) ecalRecHits = ecalRecHitsEB_;
  else      ecalRecHits = ecalRecHitsEE_;

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  Float_t rho = *rhoH;

  // Direction of the superClusters
  eta_c = superCluster.eta();
  phi_c = superCluster.phi();
    
  // Raw energy
  rawEnergy_c = superCluster.rawEnergy();
  preshowerEnergy_c = superCluster.preshowerEnergy();
  preshowerEnergyPlane1_c = superCluster.preshowerEnergyPlane1();
  preshowerEnergyPlane2_c = superCluster.preshowerEnergyPlane2();
    
  // Width of the cluster
  scEta_c = superCluster.position().Eta();
  scPhi_c = superCluster.position().Phi();
  etaWidth_c = superCluster.etaWidth();
  phiWidth_c = superCluster.phiWidth();
  delEtaSeed_c = superCluster.seed()->eta() - superCluster.position().Eta();
  delPhiSeed_c = reco::deltaPhi(superCluster.seed()->phi(),superCluster.position().Phi());
    
  // Fiducial information
  isEB_c = isEB;
  isEE_c = !isEB;

  rhoValue_c = rho;
    
  // Information about supercluster
  seedEnergy_c = superCluster.seed()->energy();
  numberOfClusters_c = superCluster.clusters().size();

  EcalClusterLocal ecalLocal;    

  clusterRawEnergy_c.clear();
  clusterPhi_c.clear();
  clusterEta_c.clear();
  clusterDPhiToSeed_c.clear();
  clusterDEtaToSeed_c.clear();
  iEtaSeed_c.clear();
  iPhiSeed_c.clear();
  iEtaCoordinate_c.clear();
  iPhiCoordinate_c.clear();
  cryEtaCoordinate_c.clear();
  cryPhiCoordinate_c.clear();
  iXSeed_c.clear();
  iYSeed_c.clear();
  iXCoordinate_c.clear();
  iYCoordinate_c.clear();
  cryXCoordinate_c.clear();
  cryYCoordinate_c.clear();
  iEtaMod5_c.clear();
  iPhiMod2_c.clear();
  iEtaMod20_c.clear();
  iPhiMod20_c.clear();
  r9_c.clear();
  sigmaIetaIeta_c.clear();
  sigmaIetaIphi_c.clear();
  sigmaIphiIphi_c.clear();
  e5x5_c.clear();
  e3x3_c.clear();
  eMax_c.clear();
  e2nd_c.clear();
  eTop_c.clear();
  eBottom_c.clear();
  eLeft_c.clear();
  eRight_c.clear();
  e2x5Max_c.clear();
  e2x5Left_c.clear();
  e2x5Right_c.clear();
  e2x5Top_c.clear();
  e2x5Bottom_c.clear();
  full5x5_r9_c.clear();
  full5x5_sigmaIetaIeta_c.clear();
  full5x5_sigmaIetaIphi_c.clear();
  full5x5_sigmaIphiIphi_c.clear();
  full5x5_e5x5_c.clear();
  full5x5_e3x3_c.clear();
  full5x5_eMax_c.clear();
  full5x5_e2nd_c.clear();
  full5x5_eTop_c.clear();
  full5x5_eBottom_c.clear();
  full5x5_eLeft_c.clear();
  full5x5_eRight_c.clear();
  full5x5_e2x5Max_c.clear();
  full5x5_e2x5Left_c.clear();
  full5x5_e2x5Right_c.clear();
  full5x5_e2x5Top_c.clear();
  full5x5_e2x5Bottom_c.clear();

  auto& cluster = *superCluster.seed();

  // Energy
  clusterRawEnergy_c.push_back(cluster.energy());
  clusterPhi_c.push_back(cluster.phi());
  clusterEta_c.push_back(cluster.eta());
  clusterDPhiToSeed_c.push_back(reco::deltaPhi(cluster.phi(), superCluster.seed()->phi()));
  clusterDEtaToSeed_c.push_back(cluster.eta() - superCluster.seed()->eta());

  if (isEB_c) {
    EBDetId ebseedid(cluster.seed());
    // Position EB
    iEtaSeed_c.push_back(ebseedid.ieta());
    iPhiSeed_c.push_back(ebseedid.iphi());

    // This is the seed, so fill saturation information
    isGap_c = EBDetId::isNextToBoundary(ebseedid) ? 1 : 0;
    isEBEEGap_c = EBDetId::isNextToBoundary(ebseedid) ? 1 : 0;
    isEBGap_c = EBDetId::isNextToBoundary(ebseedid) ? 1 : 0;
    isEBEtaGap_c = EBDetId::isNextToEtaBoundary(ebseedid) ? 1 : 0;
    isEBPhiGap_c = EBDetId::isNextToEtaBoundary(ebseedid) ? 1 : 0;
    isEEGap_c = 0;
    isEEDeeGap_c = 0;
    isEERingGap_c = 0;

    nrSaturatedCrysIn5x5_c = EcalClusterToolsT<false>::nrSaturatedCrysIn5x5(ebseedid, ecalRecHits.product(), topology_);
    nextToDead_c = EcalTools::isNextToDead(ebseedid, iSetup) ? 1 : 0;
      
    Int_t iEta, iPhi;
    Float_t cryEta, cryPhi, dummy;
    ecalLocal.localCoordsEB( cluster, iSetup,
			     cryEta, cryPhi, iEta, iPhi, dummy, dummy );
      
    iEtaCoordinate_c.push_back(iEta);
    iPhiCoordinate_c.push_back(iPhi);
    cryEtaCoordinate_c.push_back(cryEta);
    cryPhiCoordinate_c.push_back(cryPhi);
      
    Int_t iEtaSeed = ebseedid.ieta();
    Int_t iPhiSeed = ebseedid.iphi();
    Int_t signiEtaSeed = ebseedid.ieta() > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training
      
    iEtaMod5_c.push_back((iEtaSeed-signiEtaSeed)%5);
    iPhiMod2_c.push_back((iPhiSeed-1)%2);	
    iEtaMod20_c.push_back(abs(iEtaSeed)<=25 ? (iEtaSeed-signiEtaSeed) : (iEtaSeed-26*signiEtaSeed)%20);
    iPhiMod20_c.push_back((iPhiSeed-1)%20);
      
    iXSeed_c.push_back(-999);
    iYSeed_c.push_back(-999);
    iXCoordinate_c.push_back(-999);
    iYCoordinate_c.push_back(-999);
    cryXCoordinate_c.push_back(-999);
    cryYCoordinate_c.push_back(-999);
  } else {
    EEDetId eeseedid(cluster.seed());

    // This is the seed, so fill saturation information
    isGap_c = EEDetId::isNextToBoundary(eeseedid) ? 1 : 0;
    isEBEEGap_c = EEDetId::isNextToBoundary(eeseedid) ? 1 : 0;
    isEBGap_c = EEDetId::isNextToBoundary(eeseedid) ? 1 : 0;
    isEBEtaGap_c = 0;
    isEBPhiGap_c = 0;
    isEEGap_c = EEDetId::isNextToBoundary(eeseedid) ? 1 : 0;
    isEEDeeGap_c = EEDetId::isNextToDBoundary(eeseedid) ? 1 : 0;
    isEERingGap_c = EEDetId::isNextToRingBoundary(eeseedid) ? 1 : 0;

    nrSaturatedCrysIn5x5_c = EcalClusterToolsT<false>::nrSaturatedCrysIn5x5(eeseedid, ecalRecHits.product(), topology_);
    nextToDead_c = EcalTools::isNextToDead(eeseedid, iSetup) ? 1 : 0;

    // Position EE
    iXSeed_c.push_back(eeseedid.ix());
    iYSeed_c.push_back(eeseedid.iy());
    Int_t iX, iY;
    Float_t cryX, cryY, dummy;
    ecalLocal.localCoordsEE( cluster, iSetup,
			     cryX, cryY, iX, iY, dummy, dummy );
      
    iXCoordinate_c.push_back(iX);
    iYCoordinate_c.push_back(iY);
    cryXCoordinate_c.push_back(cryX);
    cryYCoordinate_c.push_back(cryY);
      
    iEtaSeed_c.push_back(-999);
    iPhiSeed_c.push_back(-999);
    iEtaCoordinate_c.push_back(-999);
    iPhiCoordinate_c.push_back(-999);
    cryEtaCoordinate_c.push_back(-999);
    cryPhiCoordinate_c.push_back(-999);
    iEtaMod5_c.push_back(-999);
    iPhiMod2_c.push_back(-999);
    iEtaMod20_c.push_back(-999);
    iPhiMod20_c.push_back(-999);
  }

  // Showershape
  std::vector<float> localCovariances = EcalClusterToolsT<false>::localCovariances(cluster, ecalRecHits.product(), topology_);
  r9_c.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
  sigmaIetaIeta_c.push_back(sqrt(localCovariances[0]));
  sigmaIetaIphi_c.push_back(localCovariances[1]/sqrt(localCovariances[0]*localCovariances[2]));
  sigmaIphiIphi_c.push_back(sqrt(localCovariances[2]));
  e5x5_c.push_back(EcalClusterToolsT<false>::e5x5(cluster, ecalRecHits.product(), topology_));
  e3x3_c.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_));
  eMax_c.push_back(EcalClusterToolsT<false>::eMax(cluster, ecalRecHits.product()));
  e2nd_c.push_back(EcalClusterToolsT<false>::e2nd(cluster, ecalRecHits.product()));
  eTop_c.push_back(EcalClusterToolsT<false>::eTop(cluster, ecalRecHits.product(), topology_));
  eBottom_c.push_back(EcalClusterToolsT<false>::eBottom(cluster, ecalRecHits.product(), topology_));
  eLeft_c.push_back(EcalClusterToolsT<false>::eLeft(cluster, ecalRecHits.product(), topology_));
  eRight_c.push_back(EcalClusterToolsT<false>::eRight(cluster, ecalRecHits.product(), topology_));
  e2x5Max_c.push_back(EcalClusterToolsT<false>::e2x5Max(cluster, ecalRecHits.product(), topology_));
  e2x5Top_c.push_back(EcalClusterToolsT<false>::e2x5Top(cluster, ecalRecHits.product(), topology_));
  e2x5Bottom_c.push_back(EcalClusterToolsT<false>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
  e2x5Left_c.push_back(EcalClusterToolsT<false>::e2x5Left(cluster, ecalRecHits.product(), topology_));
  e2x5Right_c.push_back(EcalClusterToolsT<false>::e2x5Right(cluster, ecalRecHits.product(), topology_));

  std::vector<float> full5x5_localCovariances = EcalClusterToolsT<true>::localCovariances(cluster, ecalRecHits.product(), topology_);
  full5x5_r9_c.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
  full5x5_sigmaIetaIeta_c.push_back(sqrt(full5x5_localCovariances[0]));
  full5x5_sigmaIetaIphi_c.push_back(full5x5_localCovariances[1]/sqrt(full5x5_localCovariances[0]*full5x5_localCovariances[2]));
  full5x5_sigmaIphiIphi_c.push_back(sqrt(full5x5_localCovariances[2]));
  full5x5_e5x5_c.push_back(EcalClusterToolsT<true>::e5x5(cluster, ecalRecHits.product(), topology_));
  full5x5_e3x3_c.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_));
  full5x5_eMax_c.push_back(EcalClusterToolsT<true>::eMax(cluster, ecalRecHits.product()));
  full5x5_e2nd_c.push_back(EcalClusterToolsT<true>::e2nd(cluster, ecalRecHits.product()));
  full5x5_eTop_c.push_back(EcalClusterToolsT<true>::eTop(cluster, ecalRecHits.product(), topology_));
  full5x5_eBottom_c.push_back(EcalClusterToolsT<true>::eBottom(cluster, ecalRecHits.product(), topology_));
  full5x5_eLeft_c.push_back(EcalClusterToolsT<true>::eLeft(cluster, ecalRecHits.product(), topology_));
  full5x5_eRight_c.push_back(EcalClusterToolsT<true>::eRight(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Max_c.push_back(EcalClusterToolsT<true>::e2x5Max(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Top_c.push_back(EcalClusterToolsT<true>::e2x5Top(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Bottom_c.push_back(EcalClusterToolsT<true>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Left_c.push_back(EcalClusterToolsT<true>::e2x5Left(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Right_c.push_back(EcalClusterToolsT<true>::e2x5Right(cluster, ecalRecHits.product(), topology_));

  for (auto cluster_it = superCluster.clustersBegin(); cluster_it != superCluster.clustersEnd(); ++cluster_it) {

    // Dereference to get the actual object
    auto&& clusterPtr = *cluster_it;
    auto&& cluster = *clusterPtr;

    // Continue if this is the seed
    if (clusterPtr == superCluster.seed()) continue;

    // Energy
    clusterRawEnergy_c.push_back(cluster.energy());
    clusterPhi_c.push_back(cluster.phi());
    clusterEta_c.push_back(cluster.eta());
    clusterDPhiToSeed_c.push_back(reco::deltaPhi(cluster.phi(), superCluster.seed()->phi()));
    clusterDEtaToSeed_c.push_back(cluster.eta() - superCluster.seed()->eta());

    if (isEB_c) {
      EBDetId ebseedid(cluster.seed());
      // Position EB
      iEtaSeed_c.push_back(ebseedid.ieta());
      iPhiSeed_c.push_back(ebseedid.iphi());

      Int_t iEta, iPhi;
      Float_t cryEta, cryPhi, dummy;
      ecalLocal.localCoordsEB( cluster, iSetup,
			       cryEta, cryPhi, iEta, iPhi, dummy, dummy );

      iEtaCoordinate_c.push_back(iEta);
      iPhiCoordinate_c.push_back(iPhi);
      cryEtaCoordinate_c.push_back(cryEta);
      cryPhiCoordinate_c.push_back(cryPhi);
	
      Int_t iEtaSeed = ebseedid.ieta();
      Int_t iPhiSeed = ebseedid.iphi();
      Int_t signiEtaSeed = ebseedid.ieta() > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training
      
      iEtaMod5_c.push_back((iEtaSeed-signiEtaSeed)%5);
      iPhiMod2_c.push_back((iPhiSeed-1)%2);	
      iEtaMod20_c.push_back(abs(iEtaSeed)<=25 ? (iEtaSeed-signiEtaSeed) : (iEtaSeed-26*signiEtaSeed)%20);
      iPhiMod20_c.push_back((iPhiSeed-1)%20);

      iXSeed_c.push_back(-999);
      iYSeed_c.push_back(-999);
      iXCoordinate_c.push_back(-999);
      iYCoordinate_c.push_back(-999);
      cryXCoordinate_c.push_back(-999);
      cryYCoordinate_c.push_back(-999);
    } else {
      EEDetId eeseedid(cluster.seed());
      // Position EE
      iXSeed_c.push_back(eeseedid.ix());
      iYSeed_c.push_back(eeseedid.iy());
      Int_t iX, iY;
      Float_t cryX, cryY, dummy;
      ecalLocal.localCoordsEE( cluster, iSetup,
			       cryX, cryY, iX, iY, dummy, dummy );
	
      iXCoordinate_c.push_back(iX);
      iYCoordinate_c.push_back(iY);
      cryXCoordinate_c.push_back(cryX);
      cryYCoordinate_c.push_back(cryY);

      iEtaSeed_c.push_back(-999);
      iPhiSeed_c.push_back(-999);
      iEtaCoordinate_c.push_back(-999);
      iPhiCoordinate_c.push_back(-999);
      cryEtaCoordinate_c.push_back(-999);
      cryPhiCoordinate_c.push_back(-999);
      iEtaMod5_c.push_back(-999);
      iPhiMod2_c.push_back(-999);
      iEtaMod20_c.push_back(-999);
      iPhiMod20_c.push_back(-999);
    }

    // Containement
    // Showershape
    std::vector<float> localCovariances = EcalClusterToolsT<false>::localCovariances(cluster, ecalRecHits.product(), topology_);
    r9_c.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
    sigmaIetaIeta_c.push_back(sqrt(localCovariances[0]));
    sigmaIetaIphi_c.push_back(localCovariances[1]/sqrt(localCovariances[0]*localCovariances[2]));
    sigmaIphiIphi_c.push_back(sqrt(localCovariances[2]));
    e5x5_c.push_back(EcalClusterToolsT<false>::e5x5(cluster, ecalRecHits.product(), topology_));
    e3x3_c.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_));
    eMax_c.push_back(EcalClusterToolsT<false>::eMax(cluster, ecalRecHits.product()));
    e2nd_c.push_back(EcalClusterToolsT<false>::e2nd(cluster, ecalRecHits.product()));
    eTop_c.push_back(EcalClusterToolsT<false>::eTop(cluster, ecalRecHits.product(), topology_));
    eBottom_c.push_back(EcalClusterToolsT<false>::eBottom(cluster, ecalRecHits.product(), topology_));
    eLeft_c.push_back(EcalClusterToolsT<false>::eLeft(cluster, ecalRecHits.product(), topology_));
    eRight_c.push_back(EcalClusterToolsT<false>::eRight(cluster, ecalRecHits.product(), topology_));
    e2x5Max_c.push_back(EcalClusterToolsT<false>::e2x5Max(cluster, ecalRecHits.product(), topology_));
    e2x5Top_c.push_back(EcalClusterToolsT<false>::e2x5Top(cluster, ecalRecHits.product(), topology_));
    e2x5Bottom_c.push_back(EcalClusterToolsT<false>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
    e2x5Left_c.push_back(EcalClusterToolsT<false>::e2x5Left(cluster, ecalRecHits.product(), topology_));
    e2x5Right_c.push_back(EcalClusterToolsT<false>::e2x5Right(cluster, ecalRecHits.product(), topology_));

    std::vector<float> full5x5_localCovariances = EcalClusterToolsT<true>::localCovariances(cluster, ecalRecHits.product(), topology_);
    full5x5_r9_c.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
    full5x5_sigmaIetaIeta_c.push_back(sqrt(full5x5_localCovariances[0]));
    full5x5_sigmaIetaIphi_c.push_back(full5x5_localCovariances[1]/sqrt(full5x5_localCovariances[0]*full5x5_localCovariances[2]));
    full5x5_sigmaIphiIphi_c.push_back(sqrt(full5x5_localCovariances[2]));
    full5x5_e5x5_c.push_back(EcalClusterToolsT<true>::e5x5(cluster, ecalRecHits.product(), topology_));
    full5x5_e3x3_c.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_));
    full5x5_eMax_c.push_back(EcalClusterToolsT<true>::eMax(cluster, ecalRecHits.product()));
    full5x5_e2nd_c.push_back(EcalClusterToolsT<true>::e2nd(cluster, ecalRecHits.product()));
    full5x5_eTop_c.push_back(EcalClusterToolsT<true>::eTop(cluster, ecalRecHits.product(), topology_));
    full5x5_eBottom_c.push_back(EcalClusterToolsT<true>::eBottom(cluster, ecalRecHits.product(), topology_));
    full5x5_eLeft_c.push_back(EcalClusterToolsT<true>::eLeft(cluster, ecalRecHits.product(), topology_));
    full5x5_eRight_c.push_back(EcalClusterToolsT<true>::eRight(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Max_c.push_back(EcalClusterToolsT<true>::e2x5Max(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Top_c.push_back(EcalClusterToolsT<true>::e2x5Top(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Bottom_c.push_back(EcalClusterToolsT<true>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Left_c.push_back(EcalClusterToolsT<true>::e2x5Left(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Right_c.push_back(EcalClusterToolsT<true>::e2x5Right(cluster, ecalRecHits.product(), topology_));
  }
	
  corrEnergy74X_c = superCluster.energy();
  corrEnergy74XError_c = 0.0;

  superClusterTree_->Fill();

}


bool SimpleNtuplizer::matchSuperClusterToGenParticle(const reco::SuperCluster& superCluster) {
  
  MCTruthHelper<reco::GenParticle> helper;
  // Maximum match radius
  double match_MaxDR = 0.3;

  // Keep track of minimum dX's
  double minDr   = 1e6;
  double minDe   = 1e6;
  double minDeDr = 1e6;

  // dX's of the match between current superCluster and genParticle
  double this_dr;
  double this_de;
  double this_dedr;

  // Only use the superCluster if it's matched successfully
  bool successful_match = false;
  const reco::GenParticle* matched_genParticle;

  // =====================================
  // Loop over genParticles

  for (const reco::GenParticle &genParticle : *genParticles_) {

    // Continue if pdgId is not 11 or status is not 1
    if( !(abs(genParticle.pdgId()) == 11 && genParticle.status()==1) 
	&&
	!(abs(genParticle.pdgId()) == 22 && genParticle.status()==1) 
	) continue;

    // Calculate distance variables
    this_dr   = reco::deltaR( genParticle, superCluster );
    this_de   = fabs( genParticle.energy()- superCluster.energy() ) / genParticle.energy();
    this_dedr = sqrt( this_dr*this_dr + this_de*this_de );

    if( this_dr < match_MaxDR
	&& this_dedr < minDeDr ) {

      minDr   = this_dr;
      minDe   = this_de;
      minDeDr = this_dedr;

      successful_match = true;
      matched_genParticle = &genParticle;
	  
    }
  }

  genMatchdR_e    = -999;
  genMatchdE_e    = -999;
  genMatchdRdE_e  = -999;
  genPt_e         = -999;
  genPhi_e        = -999;
  genEta_e        = -999;
  genMass_e       = -999;
  genEnergy_e     = -999;
  genPdgId_e      = -999;
  genStatus_e     = -999;

  // Return if particle could not be matched
  if(!successful_match) return successful_match;

  genMatchdR_e    = minDr;
  genMatchdE_e    = minDe;
  genMatchdRdE_e  = minDeDr;
  genPt_e         = matched_genParticle->pt();
  genPhi_e        = matched_genParticle->phi();
  genEta_e        = matched_genParticle->eta();
  genMass_e       = matched_genParticle->mass();
  genEnergy_e     = matched_genParticle->energy();
  genBornEnergy_e = helper.firstCopy(*matched_genParticle)->energy();
  genPdgId_e      = matched_genParticle->pdgId();
  genStatus_e     = matched_genParticle->status();

  // Return successful match value (should be true)
  return successful_match;

}
