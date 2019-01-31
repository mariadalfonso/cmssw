#include "SimpleHcalNtuplizer.h"

void SimpleHcalNtuplizer::setPFVariables(const edm::Event& iEvent, 
				     const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  
  ////clear all teh vector elements
  nClus_pf        = 0;
  clusrawE_pf     = -99;
  cluscorrE_pf    = -99;
  clusPt_pf       = -99;
  clusEta_pf      = -99;
  clusRho_pf      = -99;
  clusPhi_pf      = -99;
  clusLayer_pf    = -99;
  clusHits_pf     = -99;
  clusPS1_pf      = -99;
  clusPS2_pf      = -99;
  nvtx_pf         = -99;
  genEnergy_pf    = -99;
  genPt_pf        = -99;
  genEta_pf       = -99;
  genPhi_pf       = -99;
  genStatusFlag_pf= -99;
  ietamod20_pf      = -999;
  iphimod20_pf      = -999;
  tgtvar_pf         = -999;
  nlgtgtvar_pf      = -999;
  //  nhits_pf          = -999;
  
  edm::Handle<reco::PFClusterCollection> clustersH;
  iEvent.getByToken(pfLabel_,clustersH);
  

  /*
  edm::Handle<edm::ValueMap<reco::GenParticleRef> > clustergenH;
  iEvent.getByToken(genpfLabel_,clustergenH);

  edm::Handle<edm::ValueMap<int> > clusSize;
  iEvent.getByToken(clusSizeLabel_,clusSize);

  */

  edm::Handle<reco::VertexCollection> vertices;
  if(doVertex) {
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) nPV_ = 0;
    else nvtx_pf = vertices->size();
  }


  if (clustersH.isValid()) {
    // double size = (*clustersH).size();

    //    size_t iP(0);
    for (auto&& pfc : *clustersH) {
      
      ///raw energy
      clusrawE_pf = pfc.energy();
      
      ///corrected energy
      cluscorrE_pf = pfc.correctedEnergy();
      
      //std::cout<<"corrected cluster energy "<<cluscorrE_pf<<std::endl;
      ///pt
      clusPt_pf = pfc.pt();
      
      ///layer number
      int layerNum = 0;
      
      PFLayer::Layer layer = pfc.layer();

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
              
      clusLayer_pf = layerNum;
 
      ///position
      auto const & crep = pfc.position();
      double eta = crep.eta();
      double phi = crep.phi();


      unsigned int nhits=0;
      const std::vector< reco::PFRecHitFraction >& fracs = pfc.recHitFractions();
      for(unsigned i=0; i<fracs.size(); i++) {
	// PFRecHit is not available, print the detID
	if( !fracs[i].recHitRef().isAvailable() ) {
	  nhits++;
	  //	  cout<< pfc.printHitAndFraction(i) << endl;
	}
      }

      double depth = pfc.depth();
      //      double rho = crep.rho();

      //https://github.com/cms-sw/cmssw/blob/02d4198c0b6615287fd88e9a8ff650aea994412e/RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericPFlowPositionCalc.cc#L54
      //      auto nhits = pfc.recHitFractions().size();

      /*
      auto const recHitCollection = &(*pfc.recHitFractions()[0].recHitRef()) - pfc.recHitFractions()[0].recHitRef().key();

      for(auto i=0U; i<nhits; ++i) {
	auto const & hf = pfc.recHitFractions()[i];
	auto k = hf.recHitRef().key();
	auto p = recHitCollection+k;
	cout << " (*p).energy() = " << (*p).energy() << " float(hf.fraction())=" << float(hf.fraction()) << endl;

	//	hits[i]= {p,(*p).energy(), float(hf.fraction())}; 
      }      
      */

      clusEta_pf = eta;
      clusPhi_pf = phi;
      //      clusRho_pf = rho;
      clusHits_pf = nhits;

      cout << "------------------------------------------------------" << endl;

      //      cout << "layerNum=" << layerNum << " eta=" << eta << " phi=" << phi << " nhits=" << nhits << endl;

      matchClusterToGenParticle(pfc);

      float response=pfc.energy()/genEnergy_e;

      if(abs(eta)> 1.5 && pfc.energy()>50 && clusHits_pf>2 && response < 1.5  && response > 0  ) { 
	for(unsigned i=0; i<fracs.size(); i++) {
	  // PFRecHit is not available, print the detID
	  if( !fracs[i].recHitRef().isAvailable() ) {
	    cout << pfc.printHitAndFraction(i) << " " << endl;
	    const auto& id = pfc.hitsAndFractions()[i].first.rawId();
	    //	    HcalDetId(id).depth();
	    cout << "ieta= " << HcalDetId(id).ieta() << " iphi= " << HcalDetId(id).iphi() << " depth= " << HcalDetId(id).depth()  << " E(3Dclus)=" <<  pfc.energy()  << " genPdgId_e = " << genPdgId_e  << " scale = " << response << endl;
	  }
	}
      }

            
      nClus_pf++;

      pfTree_->Fill();       
    }//for (reco::PFClusterCollection::const_iterator pfc=(....))
    

    
  }//if (clustersH.isValid())
  else{
    cout<<"Handle now found!!!"<<endl;
  }
     
  
}//void SimpleNtuplizer::setPFVariables


bool SimpleHcalNtuplizer::matchClusterToGenParticle(const reco::PFCluster& pfc) {

  MCTruthHelper<reco::GenParticle> helper;
  // Maximum match radius
  double match_MaxDR = 0.1;

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

    /*
      pi+ 211
      KL 130, kS 310, K0 311, k+ 321
      n 2112
    */

    // Continue if pdgId is not 11 or status is not 1
    if( !((abs(genParticle.pdgId())==211 or abs(genParticle.pdgId())==130 or abs(genParticle.pdgId())==310 or abs(genParticle.pdgId())==311 or abs(genParticle.pdgId())==321 or abs(genParticle.pdgId())==2112) && genParticle.status()==1) ) continue;

    // Calculate distance variables
    this_dr   = reco::deltaR( genParticle, pfc );
    this_de   = fabs( pfc.energy()- pfc.energy() ) / genParticle.energy();
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

