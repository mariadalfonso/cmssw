#include "SimpleNtuplizer.h"

void SimpleNtuplizer::setPFVariables(const edm::Event& iEvent, 
				     const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
  ////clear all teh vector elements
  nClus_ = 0;
  clusE_.clear();
  clusPt_.clear();
  clusVx_.clear();
  clusVy_.clear();
  clusVz_.clear();
  clusEta_.clear();
  clusRho_.clear();
  clusPhi_.clear();
  clusLayer_.clear();
  clusNrecHits_.clear();
  
  
  //std::cout<<"inside FlatTreeMaker"<<std::endl;
  
  edm::Handle<reco::PFClusterCollection> clustersH;
  iEvent.getByToken(pfLabel_,clustersH);
  
  if (clustersH.isValid()) {
    double size = (*clustersH).size();
    //if(size>0) cout<<"size is "<<size<<endl;
    //cout<<clustersH.isValid()<<endl;
    
    for (reco::PFClusterCollection::const_iterator pfc=(*clustersH).begin(); pfc!=(*clustersH).end(); pfc++){
      
      ///energy
      double energy = (*pfc).energy();
      //cout<<"Energy is "<<energy<<endl;
      
      clusE_.push_back(energy);
      
      
      ///pt
      clusPt_.push_back( (*pfc).pt() );
      
      ///layer number
      int layerNum = 0;
      
      PFLayer::Layer layer = (*pfc).layer();
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
       
       clusLayer_.push_back(layerNum);
       
       ///position
       auto const & crep = (*pfc).positionREP();
       double eta = crep.eta();
       double phi = crep.phi();
       double rho = crep.rho();
       
       clusEta_.push_back(eta);
       clusPhi_.push_back(phi);
       clusRho_.push_back(rho);
       
       ///vertex
       double vx = (*pfc).vx();
       double vy = (*pfc).vy();
       double vz = (*pfc).vz();
       
       clusVx_.push_back(vx);
       clusVy_.push_back(vy);
       clusVz_.push_back(vz);
       
       ///no. of rechits
       int nhits = (*pfc).hitsAndFractions().size();
       clusNrecHits_.push_back(nhits);
       
       nClus_++;
       
    }//for (reco::PFClusterCollection::const_iterator pfc=(....))
    
    pfTree_->Fill();
    
  }//if (clustersH.isValid())
  else{
    cout<<"Handle now found!!!"<<endl;
  }
  
  
  
  
}//void SimpleNtuplizer::setPFVariables
