#include "RecoLocalCalo/HcalRecAlgos/interface/MahiFit_gpu.h" 

namespace hcal { namespace mahi {

__device__
MahiFit::MahiFit(float const* pshape) :
  fullTSSize_{19}, 
  fullTSofInterest_{8},
  pshape_{pshape},
  functor_{pshape,
           false,false,false,
		   1,0,0,10}
{}

__device__
void MahiFit::phase1Apply(const HBHEChannelInfo& channelData,
			  float& reconstructedEnergy,
			  float& reconstructedTime,
			  bool& useTriple, 
			  float& chi2) const {

//  assert(channelData.nSamples()==8||channelData.nSamples()==10);

  resetWorkspace();

  nnlsWork_.tsSize = channelData.nSamples();
  nnlsWork_.tsOffset = channelData.soi();
  nnlsWork_.fullTSOffset = fullTSofInterest_ - nnlsWork_.tsOffset;

  // 1 sigma time constraint
  if (channelData.hasTimeInfo()) nnlsWork_.dt=timeSigmaSiPM_;
  else nnlsWork_.dt=timeSigmaHPD_;


  //Average pedestal width (for covariance matrix constraint)
  float pedVal = 0.25*( channelData.tsPedestalWidth(0)*channelData.tsPedestalWidth(0)+
			channelData.tsPedestalWidth(1)*channelData.tsPedestalWidth(1)+
			channelData.tsPedestalWidth(2)*channelData.tsPedestalWidth(2)+
			channelData.tsPedestalWidth(3)*channelData.tsPedestalWidth(3) );

  nnlsWork_.pedConstraint.setConstant(nnlsWork_.tsSize, nnlsWork_.tsSize, pedVal);
  nnlsWork_.amplitudes.resize(nnlsWork_.tsSize);
  nnlsWork_.noiseTerms.resize(nnlsWork_.tsSize);

  float reconstructedVals[3] = {0.0, -9999, -9999};
//  std::array<float,3> reconstructedVals {{ 0.0, -9999, -9999 }};
  
  float tsTOT = 0, tstrig = 0; // in GeV
  for(unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS){
    float charge = channelData.tsRawCharge(iTS);
    float ped = channelData.tsPedestal(iTS);

    nnlsWork_.amplitudes.coeffRef(iTS) = charge - ped;

    //ADC granularity
    float noiseADC = (1./sqrt(12))*channelData.tsDFcPerADC(iTS);

    //Photostatistics
    float noisePhoto = 0;
    if ( (charge-ped)>channelData.tsPedestalWidth(iTS)) {
      noisePhoto = sqrt((charge-ped)*channelData.fcByPE());
    }

    //Electronic pedestal
    float pedWidth = channelData.tsPedestalWidth(iTS);

    //Total uncertainty from all sources
    nnlsWork_.noiseTerms.coeffRef(iTS) = noiseADC*noiseADC + noisePhoto*noisePhoto + pedWidth*pedWidth;

    tsTOT += (charge - ped)*channelData.tsGain(0);
    if( iTS==nnlsWork_.tsOffset ){
      tstrig += (charge - ped)*channelData.tsGain(0);
    }
  }

  if(tstrig >= ts4Thresh_ && tsTOT > 0) {

    useTriple=false;

    // only do pre-fit with 1 pulse if chiSq threshold is positive
    if (chiSqSwitch_>0) {
      doFit(reconstructedVals,1);
      if (reconstructedVals[2]>chiSqSwitch_) {
	doFit(reconstructedVals,0); //nbx=0 means use configured BXs
	useTriple=true;
      }
    }
    else {
      doFit(reconstructedVals,0);
      useTriple=true;
    }
  }
  else{
    reconstructedVals[0] = 0.; //energy
    reconstructedVals[1] = -9999.; //time
    reconstructedVals[2] = -9999.; //chi2
  }
  
  reconstructedEnergy = reconstructedVals[0]*channelData.tsGain(0);
  reconstructedTime = reconstructedVals[1];
  chi2 = reconstructedVals[2];

}

__device__
void MahiFit::doFit(float correctedOutput[3], int nbx) const {

  unsigned int bxSize=1;

  if (nbx==1) {
    nnlsWork_.bxOffset = 0;
  }
  else {
    bxSize = bxSizeConf_;
    nnlsWork_.bxOffset = static_cast<int>(nnlsWork_.tsOffset) >= bxOffsetConf_ ? bxOffsetConf_ : nnlsWork_.tsOffset;
  }

  nnlsWork_.nPulseTot = bxSize;

  if (dynamicPed_) nnlsWork_.nPulseTot++;
  nnlsWork_.bxs.setZero(nnlsWork_.nPulseTot);

  if (nbx==1) {
    nnlsWork_.bxs.coeffRef(0) = 0;
  }
  else {
    for (unsigned int iBX=0; iBX<bxSize; ++iBX) {
      nnlsWork_.bxs.coeffRef(iBX) = activeBXs_[iBX] - ((static_cast<int>(nnlsWork_.tsOffset) + activeBXs_[0]) >= 0 ? 0 : (nnlsWork_.tsOffset + activeBXs_[0]));
    }
  }

  nnlsWork_.maxoffset = nnlsWork_.bxs.coeffRef(bxSize-1);
  if (dynamicPed_) nnlsWork_.bxs[nnlsWork_.nPulseTot-1] = pedestalBX_;

  nnlsWork_.pulseMat.setZero(nnlsWork_.tsSize,nnlsWork_.nPulseTot);  
  nnlsWork_.pulseDerivMat.setZero(nnlsWork_.tsSize,nnlsWork_.nPulseTot);

  FullSampleVector pulseShapeArray;
  FullSampleVector pulseDerivArray;

  int offset=0;
  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    offset=nnlsWork_.bxs.coeff(iBX);

    pulseShapeArray.setZero(nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset);
    pulseDerivArray.setZero(nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset);
    nnlsWork_.pulseCovArray[iBX].setZero(nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset, nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset);


    if (offset==pedestalBX_) {
      nnlsWork_.pulseMat.col(iBX) = SampleVector::Ones(nnlsWork_.tsSize);
    }
    else {

      updatePulseShape(nnlsWork_.amplitudes.coeff(nnlsWork_.tsOffset + offset), 
		       pulseShapeArray,
		       pulseDerivArray,
		       nnlsWork_.pulseCovArray[iBX]);
      

      nnlsWork_.pulseMat.col(iBX) = pulseShapeArray.segment(nnlsWork_.maxoffset - offset, nnlsWork_.tsSize);
      nnlsWork_.pulseDerivMat.col(iBX) = pulseDerivArray.segment(nnlsWork_.maxoffset-offset, nnlsWork_.tsSize);
    }
  }

  float chiSq = minimize();

  bool foundintime = false;
  unsigned int ipulseintime = 0;

  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    if (nnlsWork_.bxs.coeff(iBX)==0) {
      ipulseintime = iBX;
      foundintime = true;
    }
  }

  if (foundintime) {
    correctedOutput[0] = nnlsWork_.ampVec.coeff(ipulseintime); //charge
    if (correctedOutput[0]!=0) {
	float arrivalTime = calculateArrivalTime();
	correctedOutput[1] = arrivalTime; //time
    }
    else correctedOutput[1] = -9999;//time

    correctedOutput[2] = chiSq; //chi2

  }
  
}

__device__
float MahiFit::minimize() const {

  nnlsWork_.invcovp.setZero(nnlsWork_.tsSize,nnlsWork_.nPulseTot);
  nnlsWork_.ampVec.setZero(nnlsWork_.nPulseTot);
  nnlsWork_.aTaMat.setZero(nnlsWork_.nPulseTot, nnlsWork_.nPulseTot);
  nnlsWork_.aTbVec.setZero(nnlsWork_.nPulseTot);

  float oldChiSq=9999;
  float chiSq=oldChiSq;

  for( int iter=1; iter<nMaxItersMin_ ; ++iter) {

    updateCov();

    if (nnlsWork_.nPulseTot>1) {
      nnls();
    }
    else {
      onePulseMinimize();
    }

    float newChiSq=calculateChiSq();
    float deltaChiSq = newChiSq - chiSq;

    if (newChiSq==oldChiSq && newChiSq<chiSq) {
      break;
    }
    oldChiSq=chiSq;
    chiSq = newChiSq;

    if (std::abs(deltaChiSq)<deltaChiSqThresh_) break;

  }

  return chiSq;

}

__device__
void MahiFit::updatePulseShape(float itQ, FullSampleVector &pulseShape, FullSampleVector &pulseDeriv,
			       FullSampleMatrix &pulseCov) const {
  
  float t0=meanTime_;

  if(applyTimeSlew_) {
    if(itQ<=1.0) t0+=tsDelay1GeV_;
    else t0+=0;
    // TODO: time slew has to be resolved eventually
    //hcalTimeSlewDelay_->delay(itQ,slewFlavor_);
  }


  double pulseN[MaxSVSize];
  double pulseM[MaxSVSize];
  double pulseP[MaxSVSize];

  for (unsigned int i=0; i<MaxSVSize; i++) {
      pulseN[i] = 0;
      pulseM[i] = 0;
      pulseP[i] = 0;
  }

  const double xx[4]={t0, 1.0, 0.0, 3};
  const double xxm[4]={-nnlsWork_.dt+t0, 1.0, 0.0, 3};
  const double xxp[4]={ nnlsWork_.dt+t0, 1.0, 0.0, 3};

//  (*pfunctor_)(&xx[0]);
  functor_.singlePulseShapeFunc(&xx[0]);
  functor_.getPulseShape(pulseN);

//  (*pfunctor_)(&xxm[0]);
  functor_.singlePulseShapeFunc(&xxm[0]);
  functor_.getPulseShape(pulseM);
  
//  (*pfunctor_)(&xxp[0]);
  functor_.singlePulseShapeFunc(&xxp[0]);
  functor_.getPulseShape(pulseP);

  //in the 2018+ case where the sample of interest (SOI) is in TS3, add an extra offset to align 
  //with previous SOI=TS4 case assumed by psfPtr_->getPulseShape()
  int delta = 4 - nnlsWork_.tsOffset;

  auto invDt = 0.25 / nnlsWork_.dt;

  for (unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS) {

    pulseShape.coeffRef(iTS+nnlsWork_.maxoffset) = pulseN[iTS+delta];
    pulseDeriv.coeffRef(iTS+nnlsWork_.maxoffset) = (pulseM[iTS+delta]-pulseP[iTS+delta])*invDt;

    pulseM[iTS] -= pulseN[iTS];
    pulseP[iTS] -= pulseN[iTS];
  }

  for (unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS) {
    for (unsigned int jTS=0; jTS<iTS+1; ++jTS) {
      
      float tmp = 0.5*( pulseP[iTS+delta]*pulseP[jTS+delta] +
			 pulseM[iTS+delta]*pulseM[jTS+delta] );

      pulseCov(iTS+nnlsWork_.maxoffset,jTS+nnlsWork_.maxoffset) += tmp;
      if(jTS!=iTS) pulseCov(jTS+nnlsWork_.maxoffset,iTS+nnlsWork_.maxoffset) += tmp;
    }
  }
  
}

__device__
void MahiFit::updateCov() const {

  SampleMatrix invCovMat;
  invCovMat.setZero(nnlsWork_.tsSize, nnlsWork_.tsSize);
  invCovMat = nnlsWork_.noiseTerms.asDiagonal();
  invCovMat +=nnlsWork_.pedConstraint;

  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    if (nnlsWork_.ampVec.coeff(iBX)==0) continue;
    
    int offset=nnlsWork_.bxs.coeff(iBX);

    if (offset==pedestalBX_) continue;		       
    else { 
      invCovMat += nnlsWork_.ampVec.coeff(iBX)*nnlsWork_.ampVec.coeff(iBX)
	*nnlsWork_.pulseCovArray[offset+nnlsWork_.bxOffset].block(nnlsWork_.maxoffset-offset, nnlsWork_.maxoffset-offset, nnlsWork_.tsSize, nnlsWork_.tsSize);
    }
  }
  
  nnlsWork_.covDecomp.compute(invCovMat);
}

__device__
float MahiFit::calculateArrivalTime() const {

  int itIndex=0;

  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    int offset=nnlsWork_.bxs.coeff(iBX);
    if (offset==0) itIndex=iBX;
  }

  PulseVector residuals = nnlsWork_.pulseMat*nnlsWork_.ampVec - nnlsWork_.amplitudes;
  PulseVector solution = nnlsWork_.pulseDerivMat.colPivHouseholderQr().solve(residuals);
  float t = solution.coeff(itIndex)/nnlsWork_.ampVec.coeff(itIndex);
  t = (t>timeLimit_) ?  timeLimit_ : 
    ((t<-timeLimit_) ? -timeLimit_ : t);

  return t;

}
  
__device__
void MahiFit::nnls() const {

  const unsigned int npulse = nnlsWork_.nPulseTot;
  const unsigned int nsamples = nnlsWork_.tsSize;

  PulseVector updateWork;
  updateWork.setZero(npulse);

  PulseVector ampvecpermtest;
  ampvecpermtest.setZero(npulse);

  nnlsWork_.invcovp = nnlsWork_.covDecomp.matrixL().solve(nnlsWork_.pulseMat);
  nnlsWork_.aTaMat = nnlsWork_.invcovp.transpose().lazyProduct(nnlsWork_.invcovp);
  nnlsWork_.aTbVec = nnlsWork_.invcovp.transpose().lazyProduct(nnlsWork_.covDecomp.matrixL().solve(nnlsWork_.amplitudes));
  
  int iter = 0;
  Index idxwmax = 0;
  float wmax = 0.0;
  float threshold = nnlsThresh_;

  nnlsWork_.nP=0;
  
  while (true) {    
    if (iter>0 || nnlsWork_.nP==0) {
      if ( nnlsWork_.nP==std::min(npulse, nsamples)) break;
      
      const unsigned int nActive = npulse - nnlsWork_.nP;
      updateWork = nnlsWork_.aTbVec - nnlsWork_.aTaMat*nnlsWork_.ampVec;
      
      Index idxwmaxprev = idxwmax;
      float wmaxprev = wmax;
      wmax = updateWork.tail(nActive).maxCoeff(&idxwmax);
      
      if (wmax<threshold || (idxwmax==idxwmaxprev && wmax==wmaxprev)) {
	break;
      }
      
      if (iter>=nMaxItersNNLS_) {
	break;
      }

      //unconstrain parameter
      Index idxp = nnlsWork_.nP + idxwmax;
      nnlsUnconstrainParameter(idxp);

    }

    while (true) {
      if (nnlsWork_.nP==0) break;     

      ampvecpermtest = nnlsWork_.ampVec;
      
      solveSubmatrix(nnlsWork_.aTaMat,nnlsWork_.aTbVec,ampvecpermtest,nnlsWork_.nP);

      //check solution
      bool positive = true;
      for (unsigned int i = 0; i < nnlsWork_.nP; ++i)
        positive &= (ampvecpermtest(i) > 0);
      if (positive) {
        nnlsWork_.ampVec.head(nnlsWork_.nP) = ampvecpermtest.head(nnlsWork_.nP);
        break;
      } 

      //update parameter vector
      Index minratioidx=0;
      
      // no realizable optimization here (because it autovectorizes!)
      float minratio = std::numeric_limits<float>::max();
      for (unsigned int ipulse=0; ipulse<nnlsWork_.nP; ++ipulse) {
	if (ampvecpermtest.coeff(ipulse)<=0.) {
	  const float c_ampvec = nnlsWork_.ampVec.coeff(ipulse);
	  const float ratio = c_ampvec/(c_ampvec-ampvecpermtest.coeff(ipulse));
	  if (ratio<minratio) {
	    minratio = ratio;
	    minratioidx = ipulse;
	  }
	}
      }
      nnlsWork_.ampVec.head(nnlsWork_.nP) += minratio*(ampvecpermtest.head(nnlsWork_.nP) - nnlsWork_.ampVec.head(nnlsWork_.nP));
      
      //avoid numerical problems with later ==0. check
      nnlsWork_.ampVec.coeffRef(minratioidx) = 0.;
      
      nnlsConstrainParameter(minratioidx);
    }
   
    ++iter;

    //adaptive convergence threshold to avoid infinite loops but still
    //ensure best value is used
    if (iter%10==0) {
      threshold *= 10.;
    }
    
  }

  
}

__device__
void MahiFit::onePulseMinimize() const {

  nnlsWork_.invcovp = nnlsWork_.covDecomp.matrixL().solve(nnlsWork_.pulseMat);

  SingleMatrix aTamatval = nnlsWork_.invcovp.transpose()*nnlsWork_.invcovp;
  SingleVector aTbvecval = nnlsWork_.invcovp.transpose()*nnlsWork_.covDecomp.matrixL().solve(nnlsWork_.amplitudes);

  nnlsWork_.ampVec.coeffRef(0) = std::max(0.f, aTbvecval.coeff(0)/aTamatval.coeff(0));


}

__device__
float MahiFit::calculateChiSq() const {
  
  return (nnlsWork_.covDecomp.matrixL().solve(nnlsWork_.pulseMat*nnlsWork_.ampVec - nnlsWork_.amplitudes)).squaredNorm();
}

/*
__device__
void MahiFit::setPulseShapeTemplate(float const* pshape) {
//void MahiFit::setPulseShapeTemplate(const HcalPulseShapes::Shape& ps,const HcalTimeSlew* hcalTimeSlewDelay) {

  pshape_ = pshape;
  functor_.assign(pshape_,
                  false,false,false,
				  1,0,0,10);
  if (!(&ps == currentPulseShape_ ))
    {

      hcalTimeSlewDelay_ = hcalTimeSlewDelay;
      tsDelay1GeV_= hcalTimeSlewDelay->delay(1.0, slewFlavor_);

      resetPulseShapeTemplate(ps);
      currentPulseShape_ = &ps;
    }
}
*/

/*
__device__
void MahiFit::resetPulseShapeTemplate(const HcalPulseShapes::Shape& ps) { 
  ++ cntsetPulseShape_;

  // only the pulse shape itself from PulseShapeFunctor is used for Mahi
  // the uncertainty terms calculated inside PulseShapeFunctor are used for Method 2 only
  psfPtr_.reset(new FitterFuncs::PulseShapeFunctor(ps,false,false,false,
						   1,0,0,10));

}
*/

__device__
void MahiFit::nnlsUnconstrainParameter(Index idxp) const {
  nnlsWork_.aTaMat.col(nnlsWork_.nP).swap(nnlsWork_.aTaMat.col(idxp));
  nnlsWork_.aTaMat.row(nnlsWork_.nP).swap(nnlsWork_.aTaMat.row(idxp));
  nnlsWork_.pulseMat.col(nnlsWork_.nP).swap(nnlsWork_.pulseMat.col(idxp));
  nnlsWork_.pulseDerivMat.col(nnlsWork_.nP).swap(nnlsWork_.pulseDerivMat.col(idxp));
  Eigen::numext::swap(nnlsWork_.aTbVec.coeffRef(nnlsWork_.nP),nnlsWork_.aTbVec.coeffRef(idxp));
  Eigen::numext::swap(nnlsWork_.ampVec.coeffRef(nnlsWork_.nP),nnlsWork_.ampVec.coeffRef(idxp));
  Eigen::numext::swap(nnlsWork_.bxs.coeffRef(nnlsWork_.nP),nnlsWork_.bxs.coeffRef(idxp));
  ++nnlsWork_.nP;
}

__device__
void MahiFit::nnlsConstrainParameter(Index minratioidx) const {
  nnlsWork_.aTaMat.col(nnlsWork_.nP-1).swap(nnlsWork_.aTaMat.col(minratioidx));
  nnlsWork_.aTaMat.row(nnlsWork_.nP-1).swap(nnlsWork_.aTaMat.row(minratioidx));
  nnlsWork_.pulseMat.col(nnlsWork_.nP-1).swap(nnlsWork_.pulseMat.col(minratioidx));
  nnlsWork_.pulseDerivMat.col(nnlsWork_.nP-1).swap(nnlsWork_.pulseDerivMat.col(minratioidx));
  Eigen::numext::swap(nnlsWork_.aTbVec.coeffRef(nnlsWork_.nP-1),nnlsWork_.aTbVec.coeffRef(minratioidx));
  Eigen::numext::swap(nnlsWork_.ampVec.coeffRef(nnlsWork_.nP-1),nnlsWork_.ampVec.coeffRef(minratioidx));
  Eigen::numext::swap(nnlsWork_.bxs.coeffRef(nnlsWork_.nP-1),nnlsWork_.bxs.coeffRef(minratioidx));
  --nnlsWork_.nP;

}

__device__
void MahiFit::phase1Debug(const HBHEChannelInfo& channelData,
			  MahiDebugInfo& mdi) const {

  float recoEnergy, recoTime, chi2;
  bool use3;
  phase1Apply(channelData, recoEnergy, recoTime, use3, chi2);


  mdi.nSamples    = channelData.nSamples();
  mdi.soi         = channelData.soi();

  mdi.use3        = use3;

  mdi.inTimeConst = nnlsWork_.dt;
  mdi.inPedAvg    = 0.25*( channelData.tsPedestalWidth(0)*channelData.tsPedestalWidth(0)+
			   channelData.tsPedestalWidth(1)*channelData.tsPedestalWidth(1)+
			   channelData.tsPedestalWidth(2)*channelData.tsPedestalWidth(2)+
			   channelData.tsPedestalWidth(3)*channelData.tsPedestalWidth(3) );
  mdi.inGain      = channelData.tsGain(0);

  for (unsigned int iTS=0; iTS<channelData.nSamples(); ++iTS) {

    double charge = channelData.tsRawCharge(iTS);
    double ped = channelData.tsPedestal(iTS);

    mdi.inNoiseADC[iTS]  = (1./sqrt(12))*channelData.tsDFcPerADC(iTS);

    if ( (charge-ped)>channelData.tsPedestalWidth(iTS)) {
      mdi.inNoisePhoto[iTS] = sqrt((charge-ped)*channelData.fcByPE());
    }
    else { mdi.inNoisePhoto[iTS] = 0; }

    mdi.inPedestal[iTS]  = channelData.tsPedestalWidth(iTS);    
    mdi.totalUCNoise[iTS] = nnlsWork_.noiseTerms.coeffRef(iTS);

    if (channelData.hasTimeInfo()) {
      mdi.inputTDC[iTS] = channelData.tsRiseTime(iTS);
    }
    else { mdi.inputTDC[iTS]=-1; }

  }

  mdi.arrivalTime = recoTime;
  mdi.chiSq       = chi2;

  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    if (nnlsWork_.bxs.coeff(iBX)==0) {
      mdi.mahiEnergy=nnlsWork_.ampVec.coeff(iBX);
      for(unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS){
	mdi.count[iTS] = iTS;
	mdi.inputTS[iTS] = nnlsWork_.amplitudes.coeff(iTS);
	mdi.itPulse[iTS] = nnlsWork_.pulseMat.col(iBX).coeff(iTS);
      }
    }
    else if (nnlsWork_.bxs.coeff(iBX)==pedestalBX_) {
      mdi.pedEnergy=nnlsWork_.ampVec.coeff(iBX);
    }
    else if (nnlsWork_.bxs.coeff(iBX)==-1) {
      mdi.pEnergy=nnlsWork_.ampVec.coeff(iBX);
      for(unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS){
        mdi.pPulse[iTS] = nnlsWork_.pulseMat.col(iBX).coeff(iTS);
      }
    }
    else if (nnlsWork_.bxs.coeff(iBX)==1) {
      mdi.nEnergy=nnlsWork_.ampVec.coeff(iBX);
      for(unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS){
	mdi.nPulse[iTS] = nnlsWork_.pulseMat.col(iBX).coeff(iTS);
      }
    }
  }  
}


__device__
void MahiFit::solveSubmatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned nP) const {
  using namespace Eigen;
  switch( nP ) { // pulse matrix is always square.
  case 10:
    {
      Matrix<float,10,10> temp = mat;
      outvec.head<10>() = temp.ldlt().solve(invec.head<10>());
    }
    break;
  case 9:
    {
      Matrix<float,9,9> temp = mat.topLeftCorner<9,9>();
      outvec.head<9>() = temp.ldlt().solve(invec.head<9>());
    }
    break;
  case 8:
    {
      Matrix<float,8,8> temp = mat.topLeftCorner<8,8>();
      outvec.head<8>() = temp.ldlt().solve(invec.head<8>());
    }
    break;
  case 7:
    {
      Matrix<float,7,7> temp = mat.topLeftCorner<7,7>();
      outvec.head<7>() = temp.ldlt().solve(invec.head<7>());
    }
    break;
  case 6:
    {
      Matrix<float,6,6> temp = mat.topLeftCorner<6,6>();
      outvec.head<6>() = temp.ldlt().solve(invec.head<6>());
    }
    break;
  case 5:
    {
      Matrix<float,5,5> temp = mat.topLeftCorner<5,5>();
      outvec.head<5>() = temp.ldlt().solve(invec.head<5>());
    }
    break;
  case 4:
    {
      Matrix<float,4,4> temp = mat.topLeftCorner<4,4>();
      outvec.head<4>() = temp.ldlt().solve(invec.head<4>());
    }
    break;
  case 3: 
    {
      Matrix<float,3,3> temp = mat.topLeftCorner<3,3>();
      outvec.head<3>() = temp.ldlt().solve(invec.head<3>());
    }
    break;
  case 2:
    {
      Matrix<float,2,2> temp = mat.topLeftCorner<2,2>();
      outvec.head<2>() = temp.ldlt().solve(invec.head<2>());
    }
    break;
  case 1:
    {
      Matrix<float,1,1> temp = mat.topLeftCorner<1,1>();
      outvec.head<1>() = temp.ldlt().solve(invec.head<1>());
    }
    break;
  default:
    return;
  }
}

__device__
void MahiFit::resetWorkspace() const {

  nnlsWork_.nPulseTot=0;
  nnlsWork_.tsSize=0;
  nnlsWork_.tsOffset=0;
  nnlsWork_.fullTSOffset=0;
  nnlsWork_.bxOffset=0;
  nnlsWork_.maxoffset=0;
  nnlsWork_.dt=0;

  nnlsWork_.amplitudes.setZero();
  nnlsWork_.noiseTerms.setZero();
  nnlsWork_.pedConstraint.setZero();


}


}}
