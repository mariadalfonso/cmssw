//######################################
//# Includes
//######################################

// system include files
#include <memory>
#include <vector>
#include <iostream>

// To check types - only use if needed
#include <typeinfo>

// user include files
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "PhysicsTools/HepMCCandAlgos/interface/MCTruthHelper.h"
#include "TTree.h"
#include "Math/VectorUtil.h"
#include "RecoEgamma/EgammaTools/interface/EcalClusterLocal.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EcalDigi/interface/EBSrFlag.h"
#include "DataFormats/EcalDigi/interface/EESrFlag.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "Geometry/CaloTopology/interface/EcalTrigTowerConstituentsMap.h"
#include "DataFormats/EcalDetId/interface/EcalScDetId.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include <algorithm>
#include <vector>


////to set the status of gen particle in PF cluster code in a 16 bit integer
void setbit(UShort_t& x, UShort_t bit);


//######################################
//# Class declaration
//######################################

class SimpleNtuplizer : public edm::EDAnalyzer {

 public:
  explicit SimpleNtuplizer(const edm::ParameterSet&);

  void setElectronVariables      (const reco::GsfElectron&,  const edm::Event&, const edm::EventSetup&);
  void setPhotonVariables        (const reco::Photon&,       const edm::Event&, const edm::EventSetup&);
  void setSuperClusterVariables  (const reco::SuperCluster&, const edm::Event&, const edm::EventSetup&, bool isEB);

  void setPFVariables(const edm::Event& iEvent, const edm::EventSetup& iSetup);
  EcalTrigTowerDetId readOutUnitOf(const EBDetId& xtalId) const;
  EcalScDetId readOutUnitOf(const EEDetId& xtalId) const;
  
  bool matchElectronToGenParticle       (const reco::GsfElectron&);
  bool matchPhotonToGenParticle         (const reco::Photon&);
  bool matchSuperClusterToGenParticle   (const reco::SuperCluster&);

  bool findTag(const reco::RecoCandidate& object, const edm::Event& iEvent, const edm::EventSetup& iSetup);
	
 private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  bool doElectronTree;
  bool doPhotonTree;
  bool doSuperClusterTree;
  bool saveUnmatched;
  bool doPFTree;
  bool doVertex;
  bool doTagAndProbe;
  bool isData;

  //typedef std::vector<std::pair<unsigned long,edm::Ptr<reco::PFCluster> > > EEtoPSAssociation;
  
  // =====================================
  // Setting tokens
        
  edm::EDGetTokenT<reco::VertexCollection>            vtxToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>       genParticleToken_;
  edm::EDGetTokenT<double>                            rhoToken_; 
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> >   PUInfoToken_;
  edm::EDGetTokenT<GenEventInfoProduct>               genEvtInfoToken_;

  edm::EDGetTokenT<reco::GsfElectronCollection>       electronToken_; 
  edm::EDGetTokenT<reco::PhotonCollection>            photonToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection>      superClustersEBToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection>      superClustersEEToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit>> ecalRecHitEBToken_;
  edm::EDGetTokenT<edm::SortedCollection<EcalRecHit>> ecalRecHitEEToken_;


  //required for reading SR flags
  const EcalTrigTowerConstituentsMap * triggerTowerMap_;
  const EcalElectronicsMapping* elecMap_;
  
  const CaloTopology *topology_;
  const CaloGeometry *geometry_;

  // =====================================
  // Set some collections as class variables
  // Because they are shared in methods
  edm::Handle<reco::GenParticleCollection>       genParticles_;
  edm::Handle<edm::SortedCollection<EcalRecHit>> ecalRecHitsEB_;
  edm::Handle<edm::SortedCollection<EcalRecHit>> ecalRecHitsEE_;
  edm::Handle<std::vector<PileupSummaryInfo> >   puInfoH_;
  edm::Handle<GenEventInfoProduct>               genEvtInfo_;
	       
  edm::EDGetTokenT<reco::PFClusterCollection> pfLabel_;
  edm::EDGetTokenT<reco::PFCluster::EEtoPSAssociation> pspfLabel_;
  edm::EDGetTokenT<edm::ValueMap<reco::GenParticleRef> > genpfLabel_;
  edm::EDGetTokenT<edm::ValueMap<int> > clusSizeLabel_;

  edm::EDGetTokenT<edm::ValueMap<float> > ps1Label_;
  edm::EDGetTokenT<edm::ValueMap<float> > ps2Label_;

  edm::EDGetTokenT<EBSrFlagCollection> ebSrFlagToken_; 
  edm::EDGetTokenT<EESrFlagCollection> eeSrFlagToken_; 

  // =====================================
  // Event variables

  TTree *eventTree_;
  TTree *electronTree_;
  TTree *superClusterTree_;
  TTree *photonTree_;
  TTree *pfTree_;

  // Central event counter (specific to this output tree)
  Int_t NtupID_;

  Int_t eventNumber_;
  Int_t luminosityBlock_;
  Int_t run_;
  Float_t weight_;
  Float_t trueNumInteractions_;


  
  Int_t nPV_;
  Int_t nElectrons_;
  Int_t nElectronsMatched_;
  Int_t nPhotons_;
  Int_t nPhotonsMatched_;
  Int_t nClusters_;
  Int_t nClustersMatched_;

  // =====================================
  // Electron tree variables

  // -----------------------------
  // Variables used in training
  Int_t eventNumber_e;
  Int_t luminosityBlock_e;
  Int_t run_e;
  Int_t isMatched_e;

  // Direction of the electrons
  Float_t eta_e;
  Float_t phi_e;
  Float_t pt_e;

  // Raw energy
  Float_t rawEnergy_e;
  Float_t preshowerEnergy_e;
  Float_t preshowerEnergyPlane1_e;
  Float_t preshowerEnergyPlane2_e;

  // Position
  Float_t scEta_e;
  Float_t scPhi_e;
  Float_t etaWidth_e;
  Float_t phiWidth_e;
  Float_t delEtaSeed_e;
  Float_t delPhiSeed_e;

  // Fiducial information
  Int_t isEB_e;
  Int_t isEE_e;
  Int_t isGap_e;
  Int_t isEBEEGap_e;
  Int_t isEBGap_e;
  Int_t isEBEtaGap_e;
  Int_t isEBPhiGap_e;
  Int_t isEEGap_e;
  Int_t isEEDeeGap_e;
  Int_t isEERingGap_e;

  Float_t rhoValue_e;

  // Information about supercluster
  Float_t seedEnergy_e;
  Int_t   numberOfClusters_e;

  // Energy
  std::vector<Float_t> clusterRawEnergy_e;
  std::vector<Float_t> clusterPhi_e;
  std::vector<Float_t> clusterEta_e;
  std::vector<Float_t> clusterDPhiToSeed_e;
  std::vector<Float_t> clusterDEtaToSeed_e;
  // Position EB
  std::vector<Int_t>   iEtaSeed_e;
  std::vector<Int_t>   iPhiSeed_e;
  std::vector<Int_t>   iEtaCoordinate_e;
  std::vector<Int_t>   iPhiCoordinate_e;
  std::vector<Float_t> cryEtaCoordinate_e;
  std::vector<Float_t> cryPhiCoordinate_e;
  // Position EE
  std::vector<Int_t>   iXSeed_e;
  std::vector<Int_t>   iYSeed_e;
  std::vector<Int_t>   iXCoordinate_e;
  std::vector<Int_t>   iYCoordinate_e;
  std::vector<Float_t> cryXCoordinate_e;
  std::vector<Float_t> cryYCoordinate_e;
  // Containement
  std::vector<Int_t>   iEtaMod5_e;
  std::vector<Int_t>   iPhiMod2_e;
  std::vector<Int_t>   iEtaMod20_e;
  std::vector<Int_t>   iPhiMod20_e;
  // Showershape
  std::vector<Float_t> r9_e;
  std::vector<Float_t> sigmaIetaIeta_e;
  std::vector<Float_t> sigmaIetaIphi_e;
  std::vector<Float_t> sigmaIphiIphi_e;
  std::vector<Float_t> e5x5_e;
  std::vector<Float_t> e3x3_e;
  std::vector<Float_t> eMax_e;
  std::vector<Float_t> e2nd_e;
  std::vector<Float_t> eTop_e;
  std::vector<Float_t> eBottom_e;
  std::vector<Float_t> eLeft_e;
  std::vector<Float_t> eRight_e;
  std::vector<Float_t> e2x5Max_e;
  std::vector<Float_t> e2x5Left_e;
  std::vector<Float_t> e2x5Right_e;
  std::vector<Float_t> e2x5Top_e;
  std::vector<Float_t> e2x5Bottom_e;
  // 5x5 Showershape
  std::vector<Float_t> full5x5_r9_e;
  std::vector<Float_t> full5x5_sigmaIetaIeta_e;
  std::vector<Float_t> full5x5_sigmaIetaIphi_e;
  std::vector<Float_t> full5x5_sigmaIphiIphi_e;
  std::vector<Float_t> full5x5_e5x5_e;
  std::vector<Float_t> full5x5_e3x3_e;
  std::vector<Float_t> full5x5_eMax_e;
  std::vector<Float_t> full5x5_e2nd_e;
  std::vector<Float_t> full5x5_eTop_e;
  std::vector<Float_t> full5x5_eBottom_e;
  std::vector<Float_t> full5x5_eLeft_e;
  std::vector<Float_t> full5x5_eRight_e;
  std::vector<Float_t> full5x5_e2x5Max_e;
  std::vector<Float_t> full5x5_e2x5Left_e;
  std::vector<Float_t> full5x5_e2x5Right_e;
  std::vector<Float_t> full5x5_e2x5Top_e;
  std::vector<Float_t> full5x5_e2x5Bottom_e;
  // H/E
  
  // Now only for the seed 5x5
  Float_t hadronicOverEm_e;
  Float_t hadronic1OverEm_e;
  Float_t hadronic2OverEm_e;
  Float_t hadronicOverEmBC_e;
  Float_t hadronic1OverEmBC_e;
  Float_t hadronic2OverEmBC_e;
  Int_t nrSaturatedCrysIn5x5_e;
  Int_t nextToDead_e;

  Float_t corrEnergy74X_e;
  Float_t corrEnergy74XError_e;

  // -----------------------------
  // Ep variables (only for electron)

  Float_t trkMomentum_e;
  Float_t trkMomentumError_e;
  Float_t trkMomentumErrorEGM_e;
  Float_t trkMomentumErrorCorrected_e;
  Float_t trkMomentumRelError_e;
  Float_t trkEta_e;
  Float_t trkPhi_e;
  Float_t eOverPuncorr_e;
  Float_t ecalDriven_e;
  Float_t ecalDrivenSeed_e;
  Float_t trackerDrivenSeed_e;
  Float_t classification_e;
  Float_t fbrem_e;
  Float_t gsfchi2_e;
  Float_t gsfndof_e;
  Float_t gsfnhits_e;
	
  Float_t genMatchdR_e;
  Float_t genMatchdE_e;
  Float_t genMatchdRdE_e;
  Float_t genPt_e;
  Float_t genPhi_e;
  Float_t genEta_e;
  Float_t genMass_e;
  Float_t genEnergy_e;
  Float_t genBornEnergy_e;
  Int_t   genPdgId_e;
  Int_t   genStatus_e;

  // =====================================
  // Cluster tree variables

  // -----------------------------
  // Variables used in training
  Int_t eventNumber_c;
  Int_t luminosityBlock_c;
  Int_t run_c;
  Int_t isMatched_c;

  // Direction of the electrons
  Float_t eta_c;
  Float_t phi_c;
  Float_t pt_c;

  // Raw energy
  Float_t rawEnergy_c;
  Float_t preshowerEnergy_c;
  Float_t preshowerEnergyPlane1_c;
  Float_t preshowerEnergyPlane2_c;

  // Position
  Float_t scEta_c;
  Float_t scPhi_c;
  Float_t etaWidth_c;
  Float_t phiWidth_c;
  Float_t delEtaSeed_c;
  Float_t delPhiSeed_c;

  // Fiducial information
  Int_t isEB_c;
  Int_t isEE_c;
  Int_t isGap_c;
  Int_t isEBEEGap_c;
  Int_t isEBGap_c;
  Int_t isEBEtaGap_c;
  Int_t isEBPhiGap_c;
  Int_t isEEGap_c;
  Int_t isEEDeeGap_c;
  Int_t isEERingGap_c;

  Float_t rhoValue_c;

  // Information about supercluster
  Float_t seedEnergy_c;
  Int_t   numberOfClusters_c;

  // Energy
  std::vector<Float_t> clusterRawEnergy_c;
  std::vector<Float_t> clusterPhi_c;
  std::vector<Float_t> clusterEta_c;
  std::vector<Float_t> clusterDPhiToSeed_c;
  std::vector<Float_t> clusterDEtaToSeed_c;
  // Position EB
  std::vector<Int_t>   iEtaSeed_c;
  std::vector<Int_t>   iPhiSeed_c;
  std::vector<Int_t>   iEtaCoordinate_c;
  std::vector<Int_t>   iPhiCoordinate_c;
  std::vector<Float_t> cryEtaCoordinate_c;
  std::vector<Float_t> cryPhiCoordinate_c;
  // Position EE
  std::vector<Int_t>   iXSeed_c;
  std::vector<Int_t>   iYSeed_c;
  std::vector<Int_t>   iXCoordinate_c;
  std::vector<Int_t>   iYCoordinate_c;
  std::vector<Float_t> cryXCoordinate_c;
  std::vector<Float_t> cryYCoordinate_c;
  // Containement
  std::vector<Int_t>   iEtaMod5_c;
  std::vector<Int_t>   iPhiMod2_c;
  std::vector<Int_t>   iEtaMod20_c;
  std::vector<Int_t>   iPhiMod20_c;
  // Showershape
  std::vector<Float_t> r9_c;
  std::vector<Float_t> sigmaIetaIeta_c;
  std::vector<Float_t> sigmaIetaIphi_c;
  std::vector<Float_t> sigmaIphiIphi_c;
  std::vector<Float_t> e5x5_c;
  std::vector<Float_t> e3x3_c;
  std::vector<Float_t> eMax_c;
  std::vector<Float_t> e2nd_c;
  std::vector<Float_t> eTop_c;
  std::vector<Float_t> eBottom_c;
  std::vector<Float_t> eLeft_c;
  std::vector<Float_t> eRight_c;
  std::vector<Float_t> e2x5Max_c;
  std::vector<Float_t> e2x5Left_c;
  std::vector<Float_t> e2x5Right_c;
  std::vector<Float_t> e2x5Top_c;
  std::vector<Float_t> e2x5Bottom_c;
  // 5x5 Showershape
  std::vector<Float_t> full5x5_r9_c;
  std::vector<Float_t> full5x5_sigmaIetaIeta_c;
  std::vector<Float_t> full5x5_sigmaIetaIphi_c;
  std::vector<Float_t> full5x5_sigmaIphiIphi_c;
  std::vector<Float_t> full5x5_e5x5_c;
  std::vector<Float_t> full5x5_e3x3_c;
  std::vector<Float_t> full5x5_eMax_c;
  std::vector<Float_t> full5x5_e2nd_c;
  std::vector<Float_t> full5x5_eTop_c;
  std::vector<Float_t> full5x5_eBottom_c;
  std::vector<Float_t> full5x5_eLeft_c;
  std::vector<Float_t> full5x5_eRight_c;
  std::vector<Float_t> full5x5_e2x5Max_c;
  std::vector<Float_t> full5x5_e2x5Left_c;
  std::vector<Float_t> full5x5_e2x5Right_c;
  std::vector<Float_t> full5x5_e2x5Top_c;
  std::vector<Float_t> full5x5_e2x5Bottom_c;

  // H/E

  // Now only for the seed 5x5
  Int_t nrSaturatedCrysIn5x5_c;
  Int_t nextToDead_c;

  Float_t corrEnergy74X_c;
  Float_t corrEnergy74XError_c;

  Float_t genMatchdR_c;
  Float_t genMatchdE_c;
  Float_t genMatchdRdE_c;
  Float_t genPt_c;
  Float_t genPhi_c;
  Float_t genEta_c;
  Float_t genMass_c;
  Float_t genEnergy_c;
  Float_t genBornEnergy_c;
  Int_t   genPdgId_c;
  Int_t   genStatus_c;

  // =====================================
  // Photon tree variables

  // -----------------------------
  // Variables used in training
  Int_t eventNumber_p;
  Int_t luminosityBlock_p;
  Int_t run_p;
  Int_t isMatched_p;

  // Direction of the electrons
  Float_t eta_p;
  Float_t phi_p;
  Float_t pt_p;

  // Raw energy
  Float_t rawEnergy_p;
  Float_t preshowerEnergy_p;
  Float_t preshowerEnergyPlane1_p;
  Float_t preshowerEnergyPlane2_p;

  // Position
  Float_t scEta_p;
  Float_t scPhi_p;
  Float_t etaWidth_p;
  Float_t phiWidth_p;
  Float_t delEtaSeed_p;
  Float_t delPhiSeed_p;

  // Fiducial information
  Int_t isEB_p;
  Int_t isEE_p;
  Int_t isGap_p;
  Int_t isEBEEGap_p;
  Int_t isEBGap_p;
  Int_t isEBEtaGap_p;
  Int_t isEBPhiGap_p;
  Int_t isEEGap_p;
  Int_t isEEDeeGap_p;
  Int_t isEERingGap_p;

  Float_t rhoValue_p;

  // Information about supercluster
  Float_t seedEnergy_p;
  Int_t   numberOfClusters_p;

  // Energy
  std::vector<Float_t> clusterRawEnergy_p;
  std::vector<Float_t> clusterPhi_p;
  std::vector<Float_t> clusterEta_p;
  std::vector<Float_t> clusterDPhiToSeed_p;
  std::vector<Float_t> clusterDEtaToSeed_p;
  // Position EB
  std::vector<Int_t>   iEtaSeed_p;
  std::vector<Int_t>   iPhiSeed_p;
  std::vector<Int_t>   iEtaCoordinate_p;
  std::vector<Int_t>   iPhiCoordinate_p;
  std::vector<Float_t> cryEtaCoordinate_p;
  std::vector<Float_t> cryPhiCoordinate_p;
  // Position EE
  std::vector<Int_t>   iXSeed_p;
  std::vector<Int_t>   iYSeed_p;
  std::vector<Int_t>   iXCoordinate_p;
  std::vector<Int_t>   iYCoordinate_p;
  std::vector<Float_t> cryXCoordinate_p;
  std::vector<Float_t> cryYCoordinate_p;
  // Containement
  std::vector<Int_t>   iEtaMod5_p;
  std::vector<Int_t>   iPhiMod2_p;
  std::vector<Int_t>   iEtaMod20_p;
  std::vector<Int_t>   iPhiMod20_p;
  // Showershape
  std::vector<Float_t> r9_p;
  std::vector<Float_t> sigmaIetaIeta_p;
  std::vector<Float_t> sigmaIetaIphi_p;
  std::vector<Float_t> sigmaIphiIphi_p;
  std::vector<Float_t> e5x5_p;
  std::vector<Float_t> e3x3_p;
  std::vector<Float_t> eMax_p;
  std::vector<Float_t> e2nd_p;
  std::vector<Float_t> eTop_p;
  std::vector<Float_t> eBottom_p;
  std::vector<Float_t> eLeft_p;
  std::vector<Float_t> eRight_p;
  std::vector<Float_t> e2x5Max_p;
  std::vector<Float_t> e2x5Left_p;
  std::vector<Float_t> e2x5Right_p;
  std::vector<Float_t> e2x5Top_p;
  std::vector<Float_t> e2x5Bottom_p;
  // 5x5 Showershape
  std::vector<Float_t> full5x5_r9_p;
  std::vector<Float_t> full5x5_sigmaIetaIeta_p;
  std::vector<Float_t> full5x5_sigmaIetaIphi_p;
  std::vector<Float_t> full5x5_sigmaIphiIphi_p;
  std::vector<Float_t> full5x5_e5x5_p;
  std::vector<Float_t> full5x5_e3x3_p;
  std::vector<Float_t> full5x5_eMax_p;
  std::vector<Float_t> full5x5_e2nd_p;
  std::vector<Float_t> full5x5_eTop_p;
  std::vector<Float_t> full5x5_eBottom_p;
  std::vector<Float_t> full5x5_eLeft_p;
  std::vector<Float_t> full5x5_eRight_p;
  std::vector<Float_t> full5x5_e2x5Max_p;
  std::vector<Float_t> full5x5_e2x5Left_p;
  std::vector<Float_t> full5x5_e2x5Right_p;
  std::vector<Float_t> full5x5_e2x5Top_p;
  std::vector<Float_t> full5x5_e2x5Bottom_p;
  // H/E

  // Now only for the seed 5x5
  Float_t hadronicOverEm_p;
  Float_t hadronic1OverEm_p;
  Float_t hadronic2OverEm_p;
  Float_t hadronicOverEmBC_p;
  Float_t hadronic1OverEmBC_p;
  Float_t hadronic2OverEmBC_p;
  Int_t nrSaturatedCrysIn5x5_p;
  Int_t nextToDead_p;

  Float_t corrEnergy74X_p;
  Float_t corrEnergy74XError_p;

  Int_t isConverted_p; 

  Float_t genMatchdR_p;
  Float_t genMatchdE_p;
  Float_t genMatchdRdE_p;
  Float_t genPt_p;
  Float_t genPhi_p;
  Float_t genEta_p;
  Float_t genMass_p;
  Float_t genEnergy_p;
  Float_t genBornEnergy_p;
  Int_t   genPdgId_p;
  Int_t   genStatus_p;

  /////////////////////////
  // PF cluster tree
  // 
  Float_t    rho_pf;
  Int_t      nClus_pf;
  Float_t    clusrawE_pf;
  Float_t    cluscorrE_pf;
  Float_t    clusPt_pf;
  Float_t    clusEta_pf;
  Float_t    clusPhi_pf;
  Float_t    clusRho_pf;
  Float_t    clusLayer_pf;
  Int_t      clusSize_pf;
  Int_t      clusIetaIx_pf;
  Int_t      clusIphiIy_pf;
  Float_t    clusPS1_pf;
  Float_t    clusPS2_pf;
  Int_t      clusFlag_pf;
  Int_t      nvtx_pf;
  Float_t    genEnergy_pf;
  Float_t    genPt_pf;
  Float_t    genEta_pf;
  Float_t    genPhi_pf;
  UShort_t   genStatusFlag_pf;

  Int_t     ietamod20_pf;
  Int_t     iphimod20_pf;
  Int_t     nhits_pf;
  Float_t tgtvar_pf;
  Float_t nlgtgtvar_pf;

  //////////////////////////
  // Tag and Probe  
  //
  Float_t tp_mll;
  Float_t tp_ptll;
  Float_t tp_tagpt;
  Float_t tp_tageta;
  Float_t tp_tagphi;
		
};
