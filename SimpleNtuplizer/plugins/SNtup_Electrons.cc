#include "SimpleNtuplizer.h"

// Function that actually reads values from the AODSIM input file
void SimpleNtuplizer::setElectronVariables(
					   const reco::GsfElectron& electron,
					   const edm::Event& iEvent,
					   const edm::EventSetup& iSetup ){

  using namespace std;
  using namespace edm;
  using namespace reco;

  // Increase count of photons in event
  nElectrons_++;

  // Try to match to genParticle; quit function if electron is not matched
  isMatched_e = matchElectronToGenParticle(electron);
  if (!isMatched_e && !saveUnmatched) return;
  else nElectronsMatched_ += isMatched_e ? 1 : 0;
    

  // Get the right ecalRecHits collection (different for barrel and encap)
  edm::Handle<SortedCollection<EcalRecHit>> ecalRecHits;
  if (electron.isEB()) ecalRecHits = ecalRecHitsEB_;
  else                 ecalRecHits = ecalRecHitsEE_;

  edm::Handle< double > rhoH;
  iEvent.getByToken(rhoToken_,rhoH);
  Float_t rho = *rhoH;

  const auto& superCluster = *electron.superCluster();

  // Direction of the electrons
  eta_e = electron.eta();
  phi_e = electron.phi();
    
  // Raw energy
  rawEnergy_e = superCluster.rawEnergy();
  preshowerEnergy_e = superCluster.preshowerEnergy();
  preshowerEnergyPlane1_e = superCluster.preshowerEnergyPlane1();
  preshowerEnergyPlane2_e = superCluster.preshowerEnergyPlane2();
    
  // Width of the cluster
  scEta_e = superCluster.position().Eta();
  scPhi_e = superCluster.position().Phi();
  etaWidth_e = superCluster.etaWidth();
  phiWidth_e = superCluster.phiWidth();
  delEtaSeed_e = superCluster.seed()->eta() - superCluster.position().Eta();
  delPhiSeed_e = reco::deltaPhi(superCluster.seed()->phi(),superCluster.position().Phi());
    
  // Fiducial information
  isEB_e = electron.isEB();
  isEE_e = electron.isEE();
  isGap_e = electron.isGap();
  isEBEEGap_e = electron.isEBEEGap();
  isEBGap_e = electron.isEBGap();
  isEBEtaGap_e = electron.isEBEtaGap();
  isEBPhiGap_e = electron.isEBPhiGap();
  isEEGap_e = electron.isEEGap();
  isEEDeeGap_e = electron.isEEDeeGap();
  isEERingGap_e = electron.isEERingGap();

  rhoValue_e = rho;
    
  // Information about supercluster
  seedEnergy_e = superCluster.seed()->energy();
  numberOfClusters_e = superCluster.clusters().size();

  EcalClusterLocal ecalLocal;    

  clusterRawEnergy_e.clear();
  clusterPhi_e.clear();
  clusterEta_e.clear();
  clusterDPhiToSeed_e.clear();
  clusterDEtaToSeed_e.clear();
  iEtaSeed_e.clear();
  iPhiSeed_e.clear();
  iEtaCoordinate_e.clear();
  iPhiCoordinate_e.clear();
  cryEtaCoordinate_e.clear();
  cryPhiCoordinate_e.clear();
  iXSeed_e.clear();
  iYSeed_e.clear();
  iXCoordinate_e.clear();
  iYCoordinate_e.clear();
  cryXCoordinate_e.clear();
  cryYCoordinate_e.clear();
  iEtaMod5_e.clear();
  iPhiMod2_e.clear();
  iEtaMod20_e.clear();
  iPhiMod20_e.clear();
  r9_e.clear();
  sigmaIetaIeta_e.clear();
  sigmaIetaIphi_e.clear();
  sigmaIphiIphi_e.clear();
  e5x5_e.clear();
  e3x3_e.clear();
  eMax_e.clear();
  e2nd_e.clear();
  eTop_e.clear();
  eBottom_e.clear();
  eLeft_e.clear();
  eRight_e.clear();
  e2x5Max_e.clear();
  e2x5Left_e.clear();
  e2x5Right_e.clear();
  e2x5Top_e.clear();
  e2x5Bottom_e.clear();
  full5x5_r9_e.clear();
  full5x5_sigmaIetaIeta_e.clear();
  full5x5_sigmaIetaIphi_e.clear();
  full5x5_sigmaIphiIphi_e.clear();
  full5x5_e5x5_e.clear();
  full5x5_e3x3_e.clear();
  full5x5_eMax_e.clear();
  full5x5_e2nd_e.clear();
  full5x5_eTop_e.clear();
  full5x5_eBottom_e.clear();
  full5x5_eLeft_e.clear();
  full5x5_eRight_e.clear();
  full5x5_e2x5Max_e.clear();
  full5x5_e2x5Left_e.clear();
  full5x5_e2x5Right_e.clear();
  full5x5_e2x5Top_e.clear();
  full5x5_e2x5Bottom_e.clear();

  auto& cluster = *superCluster.seed();

  // Energy
  clusterRawEnergy_e.push_back(cluster.energy());
  clusterPhi_e.push_back(cluster.phi());
  clusterEta_e.push_back(cluster.eta());
  clusterDPhiToSeed_e.push_back(reco::deltaPhi(cluster.phi(), superCluster.seed()->phi()));
  clusterDEtaToSeed_e.push_back(cluster.eta() - superCluster.seed()->eta());

  if (isEB_e) {
    EBDetId ebseedid(cluster.seed());
    // Position EB
    iEtaSeed_e.push_back(ebseedid.ieta());
    iPhiSeed_e.push_back(ebseedid.iphi());

    // This is the seed, so fill saturation information
    nrSaturatedCrysIn5x5_e = EcalClusterToolsT<false>::nrSaturatedCrysIn5x5(ebseedid, ecalRecHits.product(), topology_);
    nextToDead_e = EcalTools::isNextToDead(ebseedid, iSetup) ? 1 : 0;
      
    Int_t iEta, iPhi;
    Float_t cryEta, cryPhi, dummy;
    ecalLocal.localCoordsEB( cluster, iSetup,
			     cryEta, cryPhi, iEta, iPhi, dummy, dummy );
      
    iEtaCoordinate_e.push_back(iEta);
    iPhiCoordinate_e.push_back(iPhi);
    cryEtaCoordinate_e.push_back(cryEta);
    cryPhiCoordinate_e.push_back(cryPhi);
      
    Int_t iEtaSeed = ebseedid.ieta();
    Int_t iPhiSeed = ebseedid.iphi();
    Int_t signiEtaSeed = ebseedid.ieta() > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training
      
    iEtaMod5_e.push_back((iEtaSeed-signiEtaSeed)%5);
    iPhiMod2_e.push_back((iPhiSeed-1)%2);	
    iEtaMod20_e.push_back(abs(iEtaSeed)<=25 ? (iEtaSeed-signiEtaSeed) : (iEtaSeed-26*signiEtaSeed)%20);
    iPhiMod20_e.push_back((iPhiSeed-1)%20);
      
    iXSeed_e.push_back(-999);
    iYSeed_e.push_back(-999);
    iXCoordinate_e.push_back(-999);
    iYCoordinate_e.push_back(-999);
    cryXCoordinate_e.push_back(-999);
    cryYCoordinate_e.push_back(-999);
  } else {
    EEDetId eeseedid(cluster.seed());

    // This is the seed, so fill saturation information
    nrSaturatedCrysIn5x5_e = EcalClusterToolsT<false>::nrSaturatedCrysIn5x5(eeseedid, ecalRecHits.product(), topology_);
    nextToDead_e = EcalTools::isNextToDead(eeseedid, iSetup) ? 1 : 0;

    // Position EE
    iXSeed_e.push_back(eeseedid.ix());
    iYSeed_e.push_back(eeseedid.iy());
    Int_t iX, iY;
    Float_t cryX, cryY, dummy;
    ecalLocal.localCoordsEE( cluster, iSetup,
			     cryX, cryY, iX, iY, dummy, dummy );
      
    iXCoordinate_e.push_back(iX);
    iYCoordinate_e.push_back(iY);
    cryXCoordinate_e.push_back(cryX);
    cryYCoordinate_e.push_back(cryY);
      
    iEtaSeed_e.push_back(-999);
    iPhiSeed_e.push_back(-999);
    iEtaCoordinate_e.push_back(-999);
    iPhiCoordinate_e.push_back(-999);
    cryEtaCoordinate_e.push_back(-999);
    cryPhiCoordinate_e.push_back(-999);
    iEtaMod5_e.push_back(-999);
    iPhiMod2_e.push_back(-999);
    iEtaMod20_e.push_back(-999);
    iPhiMod20_e.push_back(-999);
  }

  // Showershape
  std::vector<float> localCovariances = EcalClusterToolsT<false>::localCovariances(cluster, ecalRecHits.product(), topology_);
  r9_e.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
  sigmaIetaIeta_e.push_back(sqrt(localCovariances[0]));
  sigmaIetaIphi_e.push_back(localCovariances[1]/sqrt(localCovariances[0]*localCovariances[2]));
  sigmaIphiIphi_e.push_back(sqrt(localCovariances[2]));
  e5x5_e.push_back(EcalClusterToolsT<false>::e5x5(cluster, ecalRecHits.product(), topology_));
  e3x3_e.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_));
  eMax_e.push_back(EcalClusterToolsT<false>::eMax(cluster, ecalRecHits.product()));
  e2nd_e.push_back(EcalClusterToolsT<false>::e2nd(cluster, ecalRecHits.product()));
  eTop_e.push_back(EcalClusterToolsT<false>::eTop(cluster, ecalRecHits.product(), topology_));
  eBottom_e.push_back(EcalClusterToolsT<false>::eBottom(cluster, ecalRecHits.product(), topology_));
  eLeft_e.push_back(EcalClusterToolsT<false>::eLeft(cluster, ecalRecHits.product(), topology_));
  eRight_e.push_back(EcalClusterToolsT<false>::eRight(cluster, ecalRecHits.product(), topology_));
  e2x5Max_e.push_back(EcalClusterToolsT<false>::e2x5Max(cluster, ecalRecHits.product(), topology_));
  e2x5Top_e.push_back(EcalClusterToolsT<false>::e2x5Top(cluster, ecalRecHits.product(), topology_));
  e2x5Bottom_e.push_back(EcalClusterToolsT<false>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
  e2x5Left_e.push_back(EcalClusterToolsT<false>::e2x5Left(cluster, ecalRecHits.product(), topology_));
  e2x5Right_e.push_back(EcalClusterToolsT<false>::e2x5Right(cluster, ecalRecHits.product(), topology_));

  std::vector<float> full5x5_localCovariances = EcalClusterToolsT<true>::localCovariances(cluster, ecalRecHits.product(), topology_);
  full5x5_r9_e.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
  full5x5_sigmaIetaIeta_e.push_back(sqrt(full5x5_localCovariances[0]));
  full5x5_sigmaIetaIphi_e.push_back(full5x5_localCovariances[1]/sqrt(full5x5_localCovariances[0]*full5x5_localCovariances[2]));
  full5x5_sigmaIphiIphi_e.push_back(sqrt(full5x5_localCovariances[2]));
  full5x5_e5x5_e.push_back(EcalClusterToolsT<true>::e5x5(cluster, ecalRecHits.product(), topology_));
  full5x5_e3x3_e.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_));
  full5x5_eMax_e.push_back(EcalClusterToolsT<true>::eMax(cluster, ecalRecHits.product()));
  full5x5_e2nd_e.push_back(EcalClusterToolsT<true>::e2nd(cluster, ecalRecHits.product()));
  full5x5_eTop_e.push_back(EcalClusterToolsT<true>::eTop(cluster, ecalRecHits.product(), topology_));
  full5x5_eBottom_e.push_back(EcalClusterToolsT<true>::eBottom(cluster, ecalRecHits.product(), topology_));
  full5x5_eLeft_e.push_back(EcalClusterToolsT<true>::eLeft(cluster, ecalRecHits.product(), topology_));
  full5x5_eRight_e.push_back(EcalClusterToolsT<true>::eRight(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Max_e.push_back(EcalClusterToolsT<true>::e2x5Max(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Top_e.push_back(EcalClusterToolsT<true>::e2x5Top(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Bottom_e.push_back(EcalClusterToolsT<true>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Left_e.push_back(EcalClusterToolsT<true>::e2x5Left(cluster, ecalRecHits.product(), topology_));
  full5x5_e2x5Right_e.push_back(EcalClusterToolsT<true>::e2x5Right(cluster, ecalRecHits.product(), topology_));

  for (auto cluster_it = superCluster.clustersBegin(); cluster_it != superCluster.clustersEnd(); ++cluster_it) {

    // Dereference to get the actual object
    auto&& clusterPtr = *cluster_it;
    auto&& cluster = *clusterPtr;

    // Continue if this is the seed
    if (clusterPtr == superCluster.seed()) continue;

    // Energy
    clusterRawEnergy_e.push_back(cluster.energy());
    clusterPhi_e.push_back(cluster.phi());
    clusterEta_e.push_back(cluster.eta());
    clusterDPhiToSeed_e.push_back(reco::deltaPhi(cluster.phi(), superCluster.seed()->phi()));
    clusterDEtaToSeed_e.push_back(cluster.eta() - superCluster.seed()->eta());

    if (isEB_e) {
      EBDetId ebseedid(cluster.seed());
      // Position EB
      iEtaSeed_e.push_back(ebseedid.ieta());
      iPhiSeed_e.push_back(ebseedid.iphi());

      Int_t iEta, iPhi;
      Float_t cryEta, cryPhi, dummy;
      ecalLocal.localCoordsEB( cluster, iSetup,
			       cryEta, cryPhi, iEta, iPhi, dummy, dummy );

      iEtaCoordinate_e.push_back(iEta);
      iPhiCoordinate_e.push_back(iPhi);
      cryEtaCoordinate_e.push_back(cryEta);
      cryPhiCoordinate_e.push_back(cryPhi);
	
      Int_t iEtaSeed = ebseedid.ieta();
      Int_t iPhiSeed = ebseedid.iphi();
      Int_t signiEtaSeed = ebseedid.ieta() > 0 ? +1 : -1; /// this is 1*abs(ieta)/ieta in original training
      
      iEtaMod5_e.push_back((iEtaSeed-signiEtaSeed)%5);
      iPhiMod2_e.push_back((iPhiSeed-1)%2);	
      iEtaMod20_e.push_back(abs(iEtaSeed)<=25 ? (iEtaSeed-signiEtaSeed) : (iEtaSeed-26*signiEtaSeed)%20);
      iPhiMod20_e.push_back((iPhiSeed-1)%20);

      iXSeed_e.push_back(-999);
      iYSeed_e.push_back(-999);
      iXCoordinate_e.push_back(-999);
      iYCoordinate_e.push_back(-999);
      cryXCoordinate_e.push_back(-999);
      cryYCoordinate_e.push_back(-999);
    } else {
      EEDetId eeseedid(cluster.seed());
      // Position EE
      iXSeed_e.push_back(eeseedid.ix());
      iYSeed_e.push_back(eeseedid.iy());
      Int_t iX, iY;
      Float_t cryX, cryY, dummy;
      ecalLocal.localCoordsEE( cluster, iSetup,
			       cryX, cryY, iX, iY, dummy, dummy );
	
      iXCoordinate_e.push_back(iX);
      iYCoordinate_e.push_back(iY);
      cryXCoordinate_e.push_back(cryX);
      cryYCoordinate_e.push_back(cryY);

      iEtaSeed_e.push_back(-999);
      iPhiSeed_e.push_back(-999);
      iEtaCoordinate_e.push_back(-999);
      iPhiCoordinate_e.push_back(-999);
      cryEtaCoordinate_e.push_back(-999);
      cryPhiCoordinate_e.push_back(-999);
      iEtaMod5_e.push_back(-999);
      iPhiMod2_e.push_back(-999);
      iEtaMod20_e.push_back(-999);
      iPhiMod20_e.push_back(-999);
    }

    // Containement
    // Showershape
    std::vector<float> localCovariances = EcalClusterToolsT<false>::localCovariances(cluster, ecalRecHits.product(), topology_);
    r9_e.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
    sigmaIetaIeta_e.push_back(sqrt(localCovariances[0]));
    sigmaIetaIphi_e.push_back(localCovariances[1]/sqrt(localCovariances[0]*localCovariances[2]));
    sigmaIphiIphi_e.push_back(sqrt(localCovariances[2]));
    e5x5_e.push_back(EcalClusterToolsT<false>::e5x5(cluster, ecalRecHits.product(), topology_));
    e3x3_e.push_back(EcalClusterToolsT<false>::e3x3(cluster, ecalRecHits.product(), topology_));
    eMax_e.push_back(EcalClusterToolsT<false>::eMax(cluster, ecalRecHits.product()));
    e2nd_e.push_back(EcalClusterToolsT<false>::e2nd(cluster, ecalRecHits.product()));
    eTop_e.push_back(EcalClusterToolsT<false>::eTop(cluster, ecalRecHits.product(), topology_));
    eBottom_e.push_back(EcalClusterToolsT<false>::eBottom(cluster, ecalRecHits.product(), topology_));
    eLeft_e.push_back(EcalClusterToolsT<false>::eLeft(cluster, ecalRecHits.product(), topology_));
    eRight_e.push_back(EcalClusterToolsT<false>::eRight(cluster, ecalRecHits.product(), topology_));
    e2x5Max_e.push_back(EcalClusterToolsT<false>::e2x5Max(cluster, ecalRecHits.product(), topology_));
    e2x5Top_e.push_back(EcalClusterToolsT<false>::e2x5Top(cluster, ecalRecHits.product(), topology_));
    e2x5Bottom_e.push_back(EcalClusterToolsT<false>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
    e2x5Left_e.push_back(EcalClusterToolsT<false>::e2x5Left(cluster, ecalRecHits.product(), topology_));
    e2x5Right_e.push_back(EcalClusterToolsT<false>::e2x5Right(cluster, ecalRecHits.product(), topology_));

    std::vector<float> full5x5_localCovariances = EcalClusterToolsT<true>::localCovariances(cluster, ecalRecHits.product(), topology_);
    full5x5_r9_e.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_)/superCluster.rawEnergy());
    full5x5_sigmaIetaIeta_e.push_back(sqrt(full5x5_localCovariances[0]));
    full5x5_sigmaIetaIphi_e.push_back(full5x5_localCovariances[1]/sqrt(full5x5_localCovariances[0]*full5x5_localCovariances[2]));
    full5x5_sigmaIphiIphi_e.push_back(sqrt(full5x5_localCovariances[2]));
    full5x5_e5x5_e.push_back(EcalClusterToolsT<true>::e5x5(cluster, ecalRecHits.product(), topology_));
    full5x5_e3x3_e.push_back(EcalClusterToolsT<true>::e3x3(cluster, ecalRecHits.product(), topology_));
    full5x5_eMax_e.push_back(EcalClusterToolsT<true>::eMax(cluster, ecalRecHits.product()));
    full5x5_e2nd_e.push_back(EcalClusterToolsT<true>::e2nd(cluster, ecalRecHits.product()));
    full5x5_eTop_e.push_back(EcalClusterToolsT<true>::eTop(cluster, ecalRecHits.product(), topology_));
    full5x5_eBottom_e.push_back(EcalClusterToolsT<true>::eBottom(cluster, ecalRecHits.product(), topology_));
    full5x5_eLeft_e.push_back(EcalClusterToolsT<true>::eLeft(cluster, ecalRecHits.product(), topology_));
    full5x5_eRight_e.push_back(EcalClusterToolsT<true>::eRight(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Max_e.push_back(EcalClusterToolsT<true>::e2x5Max(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Top_e.push_back(EcalClusterToolsT<true>::e2x5Top(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Bottom_e.push_back(EcalClusterToolsT<true>::e2x5Bottom(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Left_e.push_back(EcalClusterToolsT<true>::e2x5Left(cluster, ecalRecHits.product(), topology_));
    full5x5_e2x5Right_e.push_back(EcalClusterToolsT<true>::e2x5Right(cluster, ecalRecHits.product(), topology_));
  }

  // Now only for the seed 5x5
  hadronicOverEm_e = electron.hadronicOverEm();
  hadronic1OverEm_e = electron.hadronicOverEm1();
  hadronic2OverEm_e = electron.hadronicOverEm2();
  hadronicOverEmBC_e = electron.hcalDepth1OverEcalBc() + electron.hcalDepth2OverEcalBc();
  hadronic1OverEmBC_e = electron.hcalDepth1OverEcalBc();
  hadronic2OverEmBC_e = electron.hcalDepth2OverEcalBc();
	
  corrEnergy74X_e = electron.corrections().correctedEcalEnergy;
  corrEnergy74XError_e = electron.corrections().correctedEcalEnergyError;

  //######################################
  //# Analyze EP
  //######################################
    
  auto&& el_track = *electron.gsfTrack();

  trkMomentum_e          = el_track.pMode();
  trkEta_e               = el_track.etaMode();
  trkPhi_e               = el_track.phiMode();
  trkMomentumRelError_e  = std::abs(el_track.qoverpModeError())*trkMomentum_e;
  trkMomentumError_e     = trkMomentumRelError_e*trkMomentum_e;

  // Other flavours of track error
  trkMomentumErrorEGM_e = electron.trackMomentumError();
  trkMomentumErrorCorrected_e = electron.p4Error(electron.corrections().candidateP4Kind);

  // E over P (uncorrected!)
  eOverPuncorr_e = (rawEnergy_e+preshowerEnergy_e)/trkMomentum_e;

  gsfchi2_e              = el_track.chi2();
  gsfndof_e              = el_track.ndof();
  gsfnhits_e             = el_track.numberOfValidHits();

  ecalDriven_e           = electron.ecalDriven();
  ecalDrivenSeed_e       = electron.ecalDrivenSeed();
  trackerDrivenSeed_e    = electron.trackerDrivenSeed();
  classification_e       = int(electron.classification());
  fbrem_e                = electron.fbrem();

  // Write class variables to the output EpTree_
  electronTree_->Fill();

}


bool SimpleNtuplizer::matchElectronToGenParticle(const reco::GsfElectron& electron) {
  
  MCTruthHelper<reco::GenParticle> helper;
  // Maximum match radius
  double match_MaxDR = 0.3;

  // Keep track of minimum dX's
  double minDr   = 1e6;
  double minDe   = 1e6;
  double minDeDr = 1e6;

  // dX's of the match between current electron and genParticle
  double this_dr;
  double this_de;
  double this_dedr;

  // Only use the electron if it's matched successfully
  bool successful_match = false;
  const reco::GenParticle* matched_genParticle;

  // =====================================
  // Loop over genParticles

  for (const reco::GenParticle &genParticle : *genParticles_) {

    // Continue if pdgId is not 11 or status is not 1
    if( !(abs(genParticle.pdgId())==11 && genParticle.status()==1) ) continue;

    // Calculate distance variables
    this_dr   = reco::deltaR( genParticle, electron );
    this_de   = fabs( genParticle.energy()- electron.energy() ) / genParticle.energy();
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
