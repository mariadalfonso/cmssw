#include "SimpleNtuplizer.h"

// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setPhotonVariables(const reco::Photon& photon,
					 const edm::Event& iEvent,
					 const edm::EventSetup& iSetup ){

  using namespace std;
  using namespace edm;
  using namespace reco;

  // Increase count of photons in event
  nPhotons_++;

  // Try to match to genParticle; quit function if photon is not matched
  isMatched_p = matchPhotonToGenParticle(photon);    
  if (!isMatched_p && !saveUnmatched) return;
  else nPhotonsMatched_ += isMatched_p ? 1 : 0;
    

  // Get the right ecalRecHits collection (different for barrel and encap)
  edm::Handle<SortedCollection<EcalRecHit>> ecalRecHits;
  if (photon.isEB()) ecalRecHits = ecalRecHitsEB_;
  else               ecalRecHits = ecalRecHitsEE_;

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  Float_t rho = *rhoH;

  const auto& superCluster = *photon.superCluster();

  // Direction of the photons
  eta_p = photon.eta();
  phi_p = photon.phi();
    
  // Raw energy
  rawEnergy_p = superCluster.rawEnergy();
  preshowerEnergy_p = superCluster.preshowerEnergy();
  preshowerEnergyPlane1_p = superCluster.preshowerEnergyPlane1();
  preshowerEnergyPlane2_p = superCluster.preshowerEnergyPlane2();
    
  // Width of the cluster
  scEta_p = superCluster.position().Eta();
  scPhi_p = superCluster.position().Phi();
  etaWidth_p = superCluster.etaWidth();
  phiWidth_p = superCluster.phiWidth();
  delEtaSeed_p = superCluster.seed()->eta() - superCluster.position().Eta();
  delPhiSeed_p = reco::deltaPhi(superCluster.seed()->phi(),superCluster.position().Phi());
    
  // Fiducial information
  isEB_p = photon.isEB();
  isEE_p = photon.isEE();
  isGap_p = (photon.isEBEEGap() || photon.isEBGap() || photon.isEEGap()) ? 1 : 0;
  isEBEEGap_p = photon.isEBEEGap();
  isEBGap_p = photon.isEBGap();
  isEBEtaGap_p = photon.isEBEtaGap();
  isEBPhiGap_p = photon.isEBPhiGap();
  isEEGap_p = photon.isEEGap();
  isEEDeeGap_p = photon.isEEDeeGap();
  isEERingGap_p = photon.isEERingGap();

  rhoValue_p = rho;
    
  // Information about supercluster
  seedEnergy_p = superCluster.seed()->energy();
  numberOfClusters_p = superCluster.clusters().size();

  EcalClusterLocal ecalLocal;    

  clusterRawEnergy_p.clear();
  clusterPhi_p.clear();
  clusterEta_p.clear();
  clusterDPhiToSeed_p.clear();
  clusterDEtaToSeed_p.clear();
  iEtaSeed_p.clear();
  iPhiSeed_p.clear();
  iEtaCoordinate_p.clear();
  iPhiCoordinate_p.clear();
  cryEtaCoordinate_p.clear();
  cryPhiCoordinate_p.clear();
  iXSeed_p.clear();
  iYSeed_p.clear();
  iXCoordinate_p.clear();
  iYCoordinate_p.clear();
  cryXCoordinate_p.clear();
  cryYCoordinate_p.clear();
  iEtaMod5_p.clear();
  iPhiMod2_p.clear();
  iEtaMod20_p.clear();
  iPhiMod20_p.clear();
  r9_p.clear();
  sigmaIetaIeta_p.clear();
  sigmaIetaIphi_p.clear();
  sigmaIphiIphi_p.clear();
  e5x5_p.clear();
  e3x3_p.clear();
  eMax_p.clear();
  e2nd_p.clear();
  eTop_p.clear();
  eBottom_p.clear();
  eLeft_p.clear();
  eRight_p.clear();
  e2x5Max_p.clear();
  e2x5Left_p.clear();
  e2x5Right_p.clear();
  e2x5Top_p.clear();
  e2x5Bottom_p.clear();
  full5x5_r9_p.clear();
  full5x5_sigmaIetaIeta_p.clear();
  full5x5_sigmaIetaIphi_p.clear();
  full5x5_sigmaIphiIphi_p.clear();
  full5x5_e5x5_p.clear();
  full5x5_e3x3_p.clear();
  full5x5_eMax_p.clear();
  full5x5_e2nd_p.clear();
  full5x5_eTop_p.clear();
  full5x5_eBottom_p.clear();
  full5x5_eLeft_p.clear();
  full5x5_eRight_p.clear();
  full5x5_e2x5Max_p.clear();
  full5x5_e2x5Left_p.clear();
  full5x5_e2x5Right_p.clear();
  full5x5_e2x5Top_p.clear();
  full5x5_e2x5Bottom_p.clear();

  auto& cluster = *superCluster.seed();

  // Energy
  clusterRawEnergy_p.push_back(cluster.energy());
  clusterPhi_p.push_back(cluster.phi());
  clusterEta_p.push_back(cluster.eta());
  clusterDPhiToSeed_p.push_back(reco::deltaPhi(cluster.phi(), superCluster.seed()->phi()));
  clusterDEtaToSeed_p.push_back(cluster.eta() - superCluster.seed()->eta());

  if (isEB_p) {
    EBDetId ebseedid(cluster.seed());
    // Position EB
    iEtaSeed_p.push_back(ebseedid.ieta());
    iPhiSeed_p.push_back(ebseedid.iphi());

    // This is the seed, so fill saturation information
    nrSaturatedCrysIn5x5_p = EcalClusterToolsT<false>::nrSaturatedCrysIn5x5(ebseedid, ecalRecHits.product(), topology_);
    nextToDead_p = EcalTools::isNextToDead(ebseedid, iSetup) ? 1 : 0;
      
    Int_t iEta, iPhi;
    Float_t cryEta, cryPhi, dummy;
    ecalLocal.localCoordsEB( cluster, iSetup,
			     cryEta, cryPhi, iEta, iPhi, dummy, dummy );
      
    iEtaCoordinate_p.push_back(iEta);
    iPhiCoordinate_p.push_back(iPhi);
    cryEtaCoordinate_p.push_back(cryEta);
    cryPhiCoordinate_p.push_back(cryPhi);
      
    Int_t iEtaSeed = ebseedid.ieta();
    Int_t iPhiSeed = ebseedid.iphi();
    Int_t signiEtaSeed = ebseedid.ieta() > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training
      
    iEtaMod5_p.push_back((iEtaSeed-signiEtaSeed)%5);
    iPhiMod2_p.push_back((iPhiSeed-1)%2);	
    iEtaMod20_p.push_back(abs(iEtaSeed)<=25 ? (iEtaSeed-signiEtaSeed) : (iEtaSeed-26*signiEtaSeed)%20);
    iPhiMod20_p.push_back((iPhiSeed-1)%20);
      
    iXSeed_p.push_back(-999);
    iYSeed_p.push_back(-999);
    iXCoordinate_p.push_back(-999);
    iYCoordinate_p.push_back(-999);
    cryXCoordinate_p.push_back(-999);
    cryYCoordinate_p.push_back(-999);
  } else {
    EEDetId eeseedid(cluster.seed());

    // This is the seed, so fill saturation information
    nrSaturatedCrysIn5x5_p = EcalClusterToolsT<false>::nrSaturatedCrysIn5x5(eeseedid, ecalRecHits.product(), topology_);
    nextToDead_p = EcalTools::isNextToDead(eeseedid, iSetup) ? 1 : 0;

    // Position EE
    iXSeed_p.push_back(eeseedid.ix());
    iYSeed_p.push_back(eeseedid.iy());
    Int_t iX, iY;
    Float_t cryX, cryY, dummy;
    ecalLocal.localCoordsEE( cluster, iSetup,
			     cryX, cryY, iX, iY, dummy, dummy );
      
    iXCoordinate_p.push_back(iX);
    iYCoordinate_p.push_back(iY);
    cryXCoordinate_p.push_back(cryX);
    cryYCoordinate_p.push_back(cryY);
      
    iEtaSeed_p.push_back(-999);
    iPhiSeed_p.push_back(-999);
    iEtaCoordinate_p.push_back(-999);
    iPhiCoordinate_p.push_back(-999);
    cryEtaCoordinate_p.push_back(-999);
    cryPhiCoordinate_p.push_back(-999);
    iEtaMod5_p.push_back(-999);
    iPhiMod2_p.push_back(-999);
    iEtaMod20_p.push_back(-999);
    iPhiMod20_p.push_back(-999);
  }

  // Showershape
  std::vector<float> localCovariances = EcalClusterToolsT<false>::localCovariances(cluster, ecalRecHits.product(), topology_);
  r9_p.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
  sigmaIetaIeta_p.push_back(sqrt(localCovariances[0]));
  sigmaIetaIphi_p.push_back(localCovariances[1]/sqrt(localCovariances[0]*localCovariances[2]));
  sigmaIphiIphi_p.push_back(sqrt(localCovariances[2]));
  e5x5_p.push_back(EcalClusterToolsT<false>::e5x5(cluster, ecalRecHits.product(), topology_));
  e3x3_p.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_));
  eMax_p.push_back(EcalClusterToolsT<false>::eMax(cluster, ecalRecHits.product()));
  e2nd_p.push_back(EcalClusterToolsT<false>::e2nd(cluster, ecalRecHits.product()));
  eTop_p.push_back(EcalClusterToolsT<false>::eTop(cluster, ecalRecHits.product(), topology_));
  eBottom_p.push_back(EcalClusterToolsT<false>::eBottom(cluster, ecalRecHits.product(), topology_));
  eLeft_p.push_back(EcalClusterToolsT<false>::eLeft(cluster, ecalRecHits.product(), topology_));
  eRight_p.push_back(EcalClusterToolsT<false>::eRight(cluster, ecalRecHits.product(), topology_));
  e2x5Max_p.push_back(EcalClusterToolsT<false>::e2x5Max(cluster, ecalRecHits.product(), topology_));
  e2x5Top_p.push_back(EcalClusterToolsT<false>::e2x5Top(cluster, ecalRecHits.product(), topology_));
  e2x5Bottom_p.push_back(EcalClusterToolsT<false>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
  e2x5Left_p.push_back(EcalClusterToolsT<false>::e2x5Left(cluster, ecalRecHits.product(), topology_));
  e2x5Right_p.push_back(EcalClusterToolsT<false>::e2x5Right(cluster, ecalRecHits.product(), topology_));

  std::vector<float> full5x5_localCovariances = EcalClusterToolsT<true>::localCovariances(cluster, ecalRecHits.product(), topology_);
  full5x5_r9_p.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
  full5x5_sigmaIetaIeta_p.push_back(sqrt(full5x5_localCovariances[0]));
  full5x5_sigmaIetaIphi_p.push_back(full5x5_localCovariances[1]/sqrt(full5x5_localCovariances[0]*full5x5_localCovariances[2]));
  full5x5_sigmaIphiIphi_p.push_back(sqrt(full5x5_localCovariances[2]));
  full5x5_e5x5_p.push_back(EcalClusterToolsT<true>::e5x5(cluster, ecalRecHits.product(), topology_));
  full5x5_e3x3_p.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_));
  full5x5_eMax_p.push_back(EcalClusterToolsT<true>::eMax(cluster, ecalRecHits.product()));
  full5x5_e2nd_p.push_back(EcalClusterToolsT<true>::e2nd(cluster, ecalRecHits.product()));
  full5x5_eTop_p.push_back(EcalClusterToolsT<true>::eTop(cluster, ecalRecHits.product(), topology_));
  full5x5_eBottom_p.push_back(EcalClusterToolsT<true>::eBottom(cluster, ecalRecHits.product(), topology_));
  full5x5_eLeft_p.push_back(EcalClusterToolsT<true>::eLeft(cluster, ecalRecHits.product(), topology_));
  full5x5_eRight_p.push_back(EcalClusterToolsT<true>::eRight(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Max_p.push_back(EcalClusterToolsT<true>::e2x5Max(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Top_p.push_back(EcalClusterToolsT<true>::e2x5Top(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Bottom_p.push_back(EcalClusterToolsT<true>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Left_p.push_back(EcalClusterToolsT<true>::e2x5Left(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Right_p.push_back(EcalClusterToolsT<true>::e2x5Right(cluster, ecalRecHits.product(), topology_));

  for (auto cluster_it = superCluster.clustersBegin(); cluster_it != superCluster.clustersEnd(); ++cluster_it) {

    // Dereference to get the actual object
    auto&& clusterPtr = *cluster_it;
    auto&& cluster = *clusterPtr;

    // Continue if this is the seed
    if (clusterPtr == superCluster.seed()) continue;

    // Energy
    clusterRawEnergy_p.push_back(cluster.energy());
    clusterPhi_p.push_back(cluster.phi());
    clusterEta_p.push_back(cluster.eta());
    clusterDPhiToSeed_p.push_back(reco::deltaPhi(cluster.phi(), superCluster.seed()->phi()));
    clusterDEtaToSeed_p.push_back(cluster.eta() - superCluster.seed()->eta());

    if (isEB_p) {
      EBDetId ebseedid(cluster.seed());
      // Position EB
      iEtaSeed_p.push_back(ebseedid.ieta());
      iPhiSeed_p.push_back(ebseedid.iphi());

      Int_t iEta, iPhi;
      Float_t cryEta, cryPhi, dummy;
      ecalLocal.localCoordsEB( cluster, iSetup,
			       cryEta, cryPhi, iEta, iPhi, dummy, dummy );

      iEtaCoordinate_p.push_back(iEta);
      iPhiCoordinate_p.push_back(iPhi);
      cryEtaCoordinate_p.push_back(cryEta);
      cryPhiCoordinate_p.push_back(cryPhi);
	
      Int_t iEtaSeed = ebseedid.ieta();
      Int_t iPhiSeed = ebseedid.iphi();
      Int_t signiEtaSeed = ebseedid.ieta() > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training
      
      iEtaMod5_p.push_back((iEtaSeed-signiEtaSeed)%5);
      iPhiMod2_p.push_back((iPhiSeed-1)%2);	
      iEtaMod20_p.push_back(abs(iEtaSeed)<=25 ? (iEtaSeed-signiEtaSeed) : (iEtaSeed-26*signiEtaSeed)%20);
      iPhiMod20_p.push_back((iPhiSeed-1)%20);

      iXSeed_p.push_back(-999);
      iYSeed_p.push_back(-999);
      iXCoordinate_p.push_back(-999);
      iYCoordinate_p.push_back(-999);
      cryXCoordinate_p.push_back(-999);
      cryYCoordinate_p.push_back(-999);
    } else {
      EEDetId eeseedid(cluster.seed());
      // Position EE
      iXSeed_p.push_back(eeseedid.ix());
      iYSeed_p.push_back(eeseedid.iy());
      Int_t iX, iY;
      Float_t cryX, cryY, dummy;
      ecalLocal.localCoordsEE( cluster, iSetup,
			       cryX, cryY, iX, iY, dummy, dummy );
	
      iXCoordinate_p.push_back(iX);
      iYCoordinate_p.push_back(iY);
      cryXCoordinate_p.push_back(cryX);
      cryYCoordinate_p.push_back(cryY);

      iEtaSeed_p.push_back(-999);
      iPhiSeed_p.push_back(-999);
      iEtaCoordinate_p.push_back(-999);
      iPhiCoordinate_p.push_back(-999);
      cryEtaCoordinate_p.push_back(-999);
      cryPhiCoordinate_p.push_back(-999);
      iEtaMod5_p.push_back(-999);
      iPhiMod2_p.push_back(-999);
      iEtaMod20_p.push_back(-999);
      iPhiMod20_p.push_back(-999);
    }

    // Containement
    // Showershape
    std::vector<float> localCovariances = EcalClusterToolsT<false>::localCovariances(cluster, ecalRecHits.product(), topology_);
    r9_p.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
    sigmaIetaIeta_p.push_back(sqrt(localCovariances[0]));
    sigmaIetaIphi_p.push_back(localCovariances[1]/sqrt(localCovariances[0]*localCovariances[2]));
    sigmaIphiIphi_p.push_back(sqrt(localCovariances[2]));
    e5x5_p.push_back(EcalClusterToolsT<false>::e5x5(cluster, ecalRecHits.product(), topology_));
    e3x3_p.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_));
    eMax_p.push_back(EcalClusterToolsT<false>::eMax(cluster, ecalRecHits.product()));
    e2nd_p.push_back(EcalClusterToolsT<false>::e2nd(cluster, ecalRecHits.product()));
    eTop_p.push_back(EcalClusterToolsT<false>::eTop(cluster, ecalRecHits.product(), topology_));
    eBottom_p.push_back(EcalClusterToolsT<false>::eBottom(cluster, ecalRecHits.product(), topology_));
    eLeft_p.push_back(EcalClusterToolsT<false>::eLeft(cluster, ecalRecHits.product(), topology_));
    eRight_p.push_back(EcalClusterToolsT<false>::eRight(cluster, ecalRecHits.product(), topology_));
    e2x5Max_p.push_back(EcalClusterToolsT<false>::e2x5Max(cluster, ecalRecHits.product(), topology_));
    e2x5Top_p.push_back(EcalClusterToolsT<false>::e2x5Top(cluster, ecalRecHits.product(), topology_));
    e2x5Bottom_p.push_back(EcalClusterToolsT<false>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
    e2x5Left_p.push_back(EcalClusterToolsT<false>::e2x5Left(cluster, ecalRecHits.product(), topology_));
    e2x5Right_p.push_back(EcalClusterToolsT<false>::e2x5Right(cluster, ecalRecHits.product(), topology_));

    std::vector<float> full5x5_localCovariances = EcalClusterToolsT<true>::localCovariances(cluster, ecalRecHits.product(), topology_);
    full5x5_r9_p.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
    full5x5_sigmaIetaIeta_p.push_back(sqrt(full5x5_localCovariances[0]));
    full5x5_sigmaIetaIphi_p.push_back(full5x5_localCovariances[1]/sqrt(full5x5_localCovariances[0]*full5x5_localCovariances[2]));
    full5x5_sigmaIphiIphi_p.push_back(sqrt(full5x5_localCovariances[2]));
    full5x5_e5x5_p.push_back(EcalClusterToolsT<true>::e5x5(cluster, ecalRecHits.product(), topology_));
    full5x5_e3x3_p.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_));
    full5x5_eMax_p.push_back(EcalClusterToolsT<true>::eMax(cluster, ecalRecHits.product()));
    full5x5_e2nd_p.push_back(EcalClusterToolsT<true>::e2nd(cluster, ecalRecHits.product()));
    full5x5_eTop_p.push_back(EcalClusterToolsT<true>::eTop(cluster, ecalRecHits.product(), topology_));
    full5x5_eBottom_p.push_back(EcalClusterToolsT<true>::eBottom(cluster, ecalRecHits.product(), topology_));
    full5x5_eLeft_p.push_back(EcalClusterToolsT<true>::eLeft(cluster, ecalRecHits.product(), topology_));
    full5x5_eRight_p.push_back(EcalClusterToolsT<true>::eRight(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Max_p.push_back(EcalClusterToolsT<true>::e2x5Max(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Top_p.push_back(EcalClusterToolsT<true>::e2x5Top(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Bottom_p.push_back(EcalClusterToolsT<true>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Left_p.push_back(EcalClusterToolsT<true>::e2x5Left(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Right_p.push_back(EcalClusterToolsT<true>::e2x5Right(cluster, ecalRecHits.product(), topology_));
  }

  // Now only for the seed 5x5
  hadronicOverEm_p = photon.hadronicOverEm();
  hadronic1OverEm_p = photon.hadronicDepth1OverEm();
  hadronic2OverEm_p = photon.hadronicDepth2OverEm();
  hadronicOverEmBC_p = photon.hadTowDepth1OverEm() + photon.hadTowDepth2OverEm();
  hadronic1OverEmBC_p = photon.hadTowDepth1OverEm();
  hadronic2OverEmBC_p = photon.hadTowDepth2OverEm();

  corrEnergy74X_p       = photon.energyCorrections().regression2Energy;
  corrEnergy74XError_p  = photon.energyCorrections().regression2EnergyError;	

  isConverted_p = photon.hasConversionTracks();
  photonTree_->Fill();

}


bool SimpleNtuplizer::matchPhotonToGenParticle(const reco::Photon& photon) {
  
  MCTruthHelper<reco::GenParticle> helper;
  // Maximum match radius
  double match_MaxDR = 0.3;

  // Keep track of minimum dX's
  double minDr   = 1e6;
  double minDe   = 1e6;
  double minDeDr = 1e6;

  // dX's of the match between current photon and genParticle
  double this_dr;
  double this_de;
  double this_dedr;

  // Only use the photon if it's matched successfully
  bool successful_match = false;
  const reco::GenParticle* matched_genParticle;

  // =====================================
  // Loop over genParticles

  for (const reco::GenParticle &genParticle : *genParticles_) {

    // Continue if pdgId is not 22 or status is not 1
    if( !(abs(genParticle.pdgId())==22 && genParticle.status()==1) ) continue;

    // Calculate distance variables
    this_dr   = reco::deltaR( genParticle, photon );
    this_de   = fabs( genParticle.energy()-photon.energy() ) / genParticle.energy();
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

  genMatchdR_p    = -999;
  genMatchdE_p    = -999;
  genMatchdRdE_p  = -999;
  genPt_p         = -999;
  genPhi_p        = -999;
  genEta_p        = -999;
  genMass_p       = -999;
  genEnergy_p     = -999;
  genPdgId_p      = -999;
  genStatus_p     = -999;

  // Return if particle could not be matched
  if(!successful_match) return successful_match;

  genMatchdR_p    = minDr;
  genMatchdE_p    = minDe;
  genMatchdRdE_p  = minDeDr;
  genPt_p         = matched_genParticle->pt();
  genPhi_p        = matched_genParticle->phi();
  genEta_p        = matched_genParticle->eta();
  genMass_p       = matched_genParticle->mass();
  genEnergy_p     = matched_genParticle->energy();
  genBornEnergy_p = helper.firstCopy(*matched_genParticle)->energy();
  genPdgId_p      = matched_genParticle->pdgId();
  genStatus_p     = matched_genParticle->status();

  // Return successful match value (should be true)
  return successful_match;

}
