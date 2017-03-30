#include "SimpleNtuplizer.h"

void SimpleNtuplizer::setPFVariables(const edm::Event& iEvent, 
				     const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
  ////clear all teh vector elements
  nClus_pf        = 0;
  clusE_pf      = -99;
  clusPt_pf     = -99;
  clusEta_pf      = -99;
  clusRho_pf      = -99;
  clusPhi_pf      = -99;
  clusLayer_pf    = -99;
  clusSize_pf = -99;
  
  
  //std::cout<<"inside FlatTreeMaker"<<std::endl;
  
  edm::Handle<reco::PFClusterCollection> clustersH;
  iEvent.getByToken(pfLabel_,clustersH);
  
  if (clustersH.isValid()) {
    double size = (*clustersH).size();
    //if(size>0) cout<<"size is "<<size<<endl;
    //cout<<clustersH.isValid()<<endl;
    
    for (reco::PFClusterCollection::const_iterator pfc=(*clustersH).begin(); pfc!=(*clustersH).end(); pfc++){
      
      ///energy
      
      clusE_pf = pfc->energy();
      
      
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
       bool iseb = pfc->layer() == PFLayer::ECAL_BARREL;

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
       


       nClus_pf++;
       
    }//for (reco::PFClusterCollection::const_iterator pfc=(....))
    
    pfTree_->Fill();
    
  }//if (clustersH.isValid())
  else{
    cout<<"Handle now found!!!"<<endl;
  }
  
  
  
  
}//void SimpleNtuplizer::setPFVariables
