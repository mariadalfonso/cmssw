#include "SimpleNtuplizer.h"


namespace {
  typedef reco::PFCluster::EEtoPSAssociation::value_type EEPSPair;
  bool sortByKey(const EEPSPair& a, const EEPSPair& b) {
    return a.first < b.first;
  } 
}

void SimpleNtuplizer::setPFVariables(const edm::Event& iEvent, 
				     const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
  ////clear all teh vector elements
  nClus_pf        = 0;
  clusrawE_pf      = -99;
  cluscorrE_pf      = -99;
  clusPt_pf     = -99;
  clusEta_pf      = -99;
  clusRho_pf      = -99;
  clusPhi_pf      = -99;
  clusLayer_pf    = -99;
  clusPS1_pf = -99;
  clusPS2_pf = -99;
  
  
  //std::cout<<"inside FlatTreeMaker"<<std::endl;
  
  edm::Handle<reco::PFClusterCollection> clustersH;
  iEvent.getByToken(pfLabel_,clustersH);
  
  edm::Handle<reco::PFCluster::EEtoPSAssociation> clusterpairH;
  iEvent.getByToken(pspfLabel_,clusterpairH);


  // Selective Readout Flags
  edm::Handle<EBSrFlagCollection> ebSrFlags;
  iEvent.getByToken(ebSrFlagToken_, ebSrFlags );
  edm::Handle<EESrFlagCollection> eeSrFlags;
  iEvent.getByToken(eeSrFlagToken_, eeSrFlags );

  ///needed for reading the SR flag
  edm::ESHandle<EcalTrigTowerConstituentsMap> hTriggerTowerMap;
  iSetup.get<IdealGeometryRecord>().get(hTriggerTowerMap);
  triggerTowerMap_ = hTriggerTowerMap.product(); 
  

  //electronics map
  edm::ESHandle< EcalElectronicsMapping > ecalmapping;
  iSetup.get< EcalMappingRcd >().get(ecalmapping);
  elecMap_ = ecalmapping.product();


  if (clustersH.isValid()) {
    double size = (*clustersH).size();
    //if(size>0) cout<<"size is "<<size<<endl;
    //cout<<clustersH.isValid()<<endl;
    
    for (reco::PFClusterCollection::const_iterator pfc=(*clustersH).begin(); pfc!=(*clustersH).end(); pfc++){
      
      ///raw energy
      clusrawE_pf = pfc->energy();
      
      ///corrected energy
      cluscorrE_pf = pfc->correctedEnergy();
      
      ///pt
      clusPt_pf = pfc->pt();
      
      ///layer number
      int layerNum = 0;
      
      PFLayer::Layer layer = pfc->layer();
      if(layer==PFLayer::PS2) layerNum = -12;
      if(layer==PFLayer::PS1) layerNum = -11;
      if(layer==PFLayer::ECAL_ENDCAP) layerNum = -2;
      if(layer==PFLayer::ECAL_BARREL) layerNum = -1;
      if(layer==PFLayer::NONE) layerNum = 0;
      if(layer==PFLayer::HCAL_BARREL1) layerNum = 1;
      if(layer==PFLayer::HCAL_BARREL2) layerNum = 2;
      
      if(layer==PFLayer::HCAL_ENDCAP) layerNum = 3;
      if(layer==PFLayer::HF_EM) layerNum = 11;
      if(layer==PFLayer::HF_HAD) layerNum = 12;
       if(layer==PFLayer::HGCAL) layerNum = 13;
       
       //cout<<"layerNum is "<<layerNum<<endl;
       
       clusLayer_pf = layerNum;
       
       ///position
       auto const & crep = pfc->positionREP();
       double eta = crep.eta();
       double phi = crep.phi();
       double rho = crep.rho();
       
       clusEta_pf = eta;
       clusPhi_pf = phi;
       clusRho_pf = rho;
       


       ///ieta, iphi
       //find seed crystal indices
       bool iseb = (pfc->layer()) == (PFLayer::ECAL_BARREL);

       if (iseb) {
	 EBDetId ebseed(pfc->seed());
	 clusIetaIx_pf = ebseed.ieta();
	 clusIphiIy_pf = ebseed.iphi();
       }
       else {
	 EEDetId eeseed(pfc->seed());
	 clusIetaIx_pf = eeseed.ix();
	 clusIphiIy_pf = eeseed.iy();      
       }
       

       ///lazy tools
       EcalClusterLazyTools lazyTool(iEvent, iSetup, ecalRecHitEBToken_, ecalRecHitEEToken_);
       clusSize_pf = lazyTool.n5x5(*pfc);
       

       ///////PS energy
       //compute preshower energies for endcap clusters
       double ePS1=0, ePS2=0;
       if(!iseb) {
	 auto ee_key_val = std::make_pair(nClus_pf,edm::Ptr<reco::PFCluster>());
            const auto clustops = std::equal_range(clusterpairH->begin(),
						   clusterpairH->end(),
						   ee_key_val,
						   sortByKey);
	    for( auto i_ps = clustops.first; i_ps != clustops.second; ++i_ps) {
	      edm::Ptr<reco::PFCluster> psclus(i_ps->second);
	      switch( psclus->layer() ) {
		
	      case PFLayer::PS1:
		ePS1 += psclus->energy();
		break;
	      case PFLayer::PS2:
		ePS2 += psclus->energy();
		break;
	      default:
		break;
	      }
	    }//for( auto i_ps = clustops.first; i_ps !=..)
	    
	    clusPS1_pf = ePS1;
	    clusPS2_pf = ePS2;
	    
       }//if(!iseb)


       ////SR flags
       if(iseb){

	 EBSrFlagCollection::const_iterator srf
	   = ebSrFlags->find(readOutUnitOf((EBDetId) pfc->seed()));
	 
	 clusFlag_pf = srf->value();
       }

       if(!iseb){
	 EESrFlagCollection::const_iterator srf
	   = eeSrFlags->find(readOutUnitOf((EEDetId)pfc->seed()));
	 
	 clusFlag_pf = srf->value();   
       }

       nClus_pf++;
       
       
    }//for (reco::PFClusterCollection::const_iterator pfc=(....))
    
    pfTree_->Fill();
    
  }//if (clustersH.isValid())
  else{
    cout<<"Handle now found!!!"<<endl;
  }
  
  
  
  
}//void SimpleNtuplizer::setPFVariables

///http://cmslxr.fnal.gov/source/Validation/EcalDigis/src/EcalSelectiveReadoutValidation.cc#0668

/*
 http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_9_0_0/doc/html/de/d38/classEcalSrFlag.html
SRF_FORCED_MASK = 0x4
SRF_FULL = 3
SRF_SUPPRESS = 0
SRF_ZS1 = 1
SRF_ZS2 = 2
*/ 


EcalTrigTowerDetId SimpleNtuplizer::readOutUnitOf(const EBDetId& xtalId) const{
  return triggerTowerMap_->towerOf(xtalId);
}


EcalScDetId SimpleNtuplizer::readOutUnitOf(const EEDetId& xtalId) const{
  const EcalElectronicsId& EcalElecId = elecMap_->getElectronicsId(xtalId);
  int iDCC= EcalElecId.dccId();
  int iDccChan = EcalElecId.towerId();
  const bool ignoreSingle = true;
  const std::vector<EcalScDetId> id = elecMap_->getEcalScDetId(iDCC, iDccChan, ignoreSingle);
  return id.size()>0?id[0]:EcalScDetId();
}
