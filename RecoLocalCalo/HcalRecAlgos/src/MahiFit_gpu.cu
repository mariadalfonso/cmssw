#define EIGEN_NO_DEBUG  // kill throws in eigen code
#include "RecoLocalCalo/HcalRecAlgos/interface/MahiFit_gpu.h" 

namespace hcal { namespace mahi {

__device__
MahiFit::MahiFit(float const* pshape) :
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
			  float& chi2,
			  float* pulseNn, float* pulseMn, float* pulsePn,
			  float* pulseShapeArray, float* pulseDerivArray, float* pulseCovArray
			  ) const {

//  assert(channelData.nSamples()==8||channelData.nSamples()==10);

  resetWorkspace();

  nnlsWork_.tsOffset = channelData.soi();
  nnlsWork_.tsSize = channelData.nSamples();

  float reconstructedVals[3] = {0.f, -9999.f, -9999.f};
//  std::array<float,3> reconstructedVals {{ 0.0, -9999, -9999 }};
  
  double tsTOT = 0, tstrig = 0; // in GeV

  for(unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS){

    auto const amplitude = channelData.tsRawCharge(iTS) - channelData.tsPedestal(iTS);
    nnlsWork_.amplitudes.coeffRef(iTS) = amplitude;

    //ADC granularity
    auto const noiseADC = norm_ * channelData.tsDFcPerADC(iTS);

    //Electronic pedestal
    auto const pedWidth = channelData.tsPedestalWidth(iTS);

    //Photostatistics
    auto const noisePhoto = (amplitude > pedWidth) ? std::sqrt(amplitude * channelData.fcByPE()) : 0.f;

    //Total uncertainty from all sources
    nnlsWork_.noiseTerms.coeffRef(iTS) = noiseADC*noiseADC + noisePhoto*noisePhoto + pedWidth*pedWidth;

    tsTOT += amplitude;
    if (iTS == nnlsWork_.tsOffset)
        tstrig += amplitude;

    }

  tsTOT *= channelData.tsGain(0);
  tstrig *= channelData.tsGain(0);


  if(tstrig >= ts4Thresh_ && tsTOT > 0) {

    //Average pedestal width (for covariance matrix constraint)
    nnlsWork_.pedVal = 0.25f*( channelData.tsPedestalWidth(0)*channelData.tsPedestalWidth(0)+
			       channelData.tsPedestalWidth(1)*channelData.tsPedestalWidth(1)+
			       channelData.tsPedestalWidth(2)*channelData.tsPedestalWidth(2)+
			       channelData.tsPedestalWidth(3)*channelData.tsPedestalWidth(3) );

    // 1 sigma time constraint
    if (channelData.hasTimeInfo()) nnlsWork_.dt=timeSigmaSiPM_;
    else nnlsWork_.dt=timeSigmaHPD_;

    nnlsWork_.amplitudes.resize(nnlsWork_.tsSize);
    nnlsWork_.noiseTerms.resize(nnlsWork_.tsSize);

    useTriple=false;

    // only do pre-fit with 1 pulse if chiSq threshold is positive
    if (chiSqSwitch_>0) {
      doFit(reconstructedVals,1, pulseNn, pulseMn, pulsePn, pulseShapeArray, pulseDerivArray, pulseCovArray);
      if (reconstructedVals[2]>chiSqSwitch_) {
	doFit(reconstructedVals,0, pulseNn, pulseMn, pulsePn, pulseShapeArray, pulseDerivArray, pulseCovArray); //nbx=0 means use configured BXs
	useTriple=true;
      }
    }
    else {
      doFit(reconstructedVals,0, pulseNn, pulseMn, pulsePn, pulseShapeArray, pulseDerivArray, pulseCovArray);
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
void MahiFit::doFit(float correctedOutput[3], int nbx, float* pulseNn, float* pulseMn, float* pulsePn, float* pulseShapeVector, float* pulseDerivVector, float* pulseCovVector) const {

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

  nnlsWork_.maxoffset = nnlsWork_.bxs.coeff(bxSize-1);
  if (dynamicPed_) nnlsWork_.bxs[nnlsWork_.nPulseTot-1] = pedestalBX_;

  //  nnlsWork_.pulseMat.setZero(nnlsWork_.tsSize,nnlsWork_.nPulseTot);
  //  if(calculateArrivalTime_) nnlsWork_.pulseDerivMat.setZero(nnlsWork_.tsSize,nnlsWork_.nPulseTot);


  const int idx = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride = blockDim.x*gridDim.x;
  if (idx >= stride ) return;

  SamplePulseMatrixMAP pulseMat_(pulseShapeVector + idx, DynStride(stride * MaxSVSize, stride));
  SamplePulseMatrixMAP pulseDerivMat_(pulseDerivVector + idx, DynStride(stride * MaxSVSize , stride));

  int sizeSQ=nnlsWork_.tsSize*nnlsWork_.tsSize;

  SampleMatrixMAP covs_[MaxPVSize] = {SampleMatrixMAP(pulseCovVector + idx, DynStride(stride * MaxSVSize, stride)),
				      SampleMatrixMAP(pulseCovVector + sizeSQ + idx, DynStride(stride * MaxSVSize, stride)),
				      SampleMatrixMAP(pulseCovVector + 2*sizeSQ +idx, DynStride(stride * MaxSVSize, stride)),
				      SampleMatrixMAP(pulseCovVector + 3*sizeSQ +idx, DynStride(stride * MaxSVSize, stride)),
				      SampleMatrixMAP(pulseCovVector + 4*sizeSQ +idx, DynStride(stride * MaxSVSize, stride)),
				      SampleMatrixMAP(pulseCovVector + 5*sizeSQ +idx, DynStride(stride * MaxSVSize, stride)),
				      SampleMatrixMAP(pulseCovVector + 6*sizeSQ +idx, DynStride(stride * MaxSVSize, stride)),
				      SampleMatrixMAP(pulseCovVector + 7*sizeSQ +idx, DynStride(stride * MaxSVSize, stride))  };

  FullSampleVector pulseShapeArray;
  FullSampleVector pulseDerivArray;
  FullSampleMatrix pulseCov;

  int offset=0;
  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    offset=nnlsWork_.bxs.coeff(iBX);

    if (offset==pedestalBX_) {
      //      nnlsWork_.pulseMat.col(iBX) = SampleVector::Ones(nnlsWork_.tsSize);
      //      if(calculateArrivalTime_) nnlsWork_.pulseDerivMat.col(iBX) = SampleVector::Zero(nnlsWork_.tsSize);
    }
    else {

      pulseShapeArray.setZero();
      if(calculateArrivalTime_) pulseDerivArray.setZero();
      pulseCov.setZero();

      //      pulseShapeArray.setZero(nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset);
      //      pulseDerivArray.setZero(nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset);
      //      pulseCov.setZero(nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset,
      //      	     		nnlsWork_.tsSize + nnlsWork_.maxoffset + nnlsWork_.bxOffset);
      //      nnlsWork_.pulseCovArray[iBX].setZero(nnlsWork_.tsSize, nnlsWork_.tsSize);

      updatePulseShape(nnlsWork_.amplitudes.coeff(nnlsWork_.tsOffset + offset), 
		       pulseShapeArray,
		       pulseDerivArray,
		       pulseCov,
		       pulseNn, pulseMn, pulsePn
		       );
      

      //      nnlsWork_.pulseMat.col(iBX) = pulseShapeArray.segment(nnlsWork_.maxoffset - offset, nnlsWork_.tsSize);
      //      if(calculateArrivalTime_) nnlsWork_.pulseDerivMat.col(iBX) = pulseDerivArray.segment(nnlsWork_.maxoffset-offset, nnlsWork_.tsSize);
      //      nnlsWork_.pulseCovArray[iBX] = pulseCov.block(
      //         			    nnlsWork_.maxoffset - offset, nnlsWork_.maxoffset - offset, nnlsWork_.tsSize, nnlsWork_.tsSize);
      pulseMat_.col(iBX) = pulseShapeArray.segment(nnlsWork_.maxoffset - offset, nnlsWork_.tsSize);
      pulseDerivMat_.col(iBX) = pulseDerivArray.segment(nnlsWork_.maxoffset - offset, nnlsWork_.tsSize);
      covs_[iBX] = pulseCov.block(nnlsWork_.maxoffset - offset, nnlsWork_.maxoffset - offset, nnlsWork_.tsSize, nnlsWork_.tsSize);

    }
  }

  const float chiSq = minimize(pulseMat_,covs_);

  bool foundintime = false;
  unsigned int ipulseintime = 0;

  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    if (nnlsWork_.bxs.coeff(iBX)==0) {
      ipulseintime = iBX;
      foundintime = true;
      break;
    }
  }

  if (foundintime) {
    correctedOutput[0] = nnlsWork_.ampVec.coeff(ipulseintime); //charge
    if (correctedOutput[0]!=0) {
        float arrivalTime = 0.f;
	if(calculateArrivalTime_) arrivalTime = calculateArrivalTime(pulseMat_,pulseDerivMat_,ipulseintime);
	correctedOutput[1] = arrivalTime; //time
    }
    else correctedOutput[1] = -9999.f;//time

    correctedOutput[2] = chiSq; //chi2

  }
}

__device__
const float MahiFit::minimize(SamplePulseMatrixMAP & pulseMat_, SampleMatrixMAP const * covs_) const {

  nnlsWork_.invcovp.setZero(nnlsWork_.tsSize,nnlsWork_.nPulseTot);
  nnlsWork_.ampVec.setZero(nnlsWork_.nPulseTot);
  //  nnlsWork_.aTaMat.setZero(nnlsWork_.nPulseTot, nnlsWork_.nPulseTot);
  //  nnlsWork_.aTbVec.setZero(nnlsWork_.nPulseTot);

  double oldChiSq=9999;
  double chiSq=oldChiSq;

  //  SampleMatrix invCovMat;
  //  invCovMat.setConstant(nnlsWork_.tsSize, nnlsWork_.tsSize, nnlsWork_.pedVal);
  //  invCovMat += nnlsWork_.noiseTerms.asDiagonal();

  for( int iter=1; iter<nMaxItersMin_ ; ++iter) {

    //    updateCov(invCovMat,covs_);
    updateCov(covs_);

    if (nnlsWork_.nPulseTot>1) {
      nnls(pulseMat_);
    }
    else {
      onePulseMinimize(pulseMat_);
    }

    const float newChiSq=calculateChiSq(pulseMat_);
    const float deltaChiSq = newChiSq - chiSq;

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
void MahiFit::updatePulseShape(const float itQ, FullSampleVector &pulseShape, FullSampleVector &pulseDeriv,
			       FullSampleMatrix &pulseCov,
			       float* pulseNn, float* pulseMn, float* pulsePn
			       ) const {
  
  const int idx = threadIdx.x + blockIdx.x * blockDim.x;
  const int stride = blockDim.x*gridDim.x;

  float t0=meanTime_;

  if(applyTimeSlew_) {
    if(itQ<=1.0) t0+=tsDelay1GeV_;
    else t0+=0.f;
    // TODO: time slew has to be resolved eventually
    //hcalTimeSlewDelay_->delay(itQ,slewFlavor_);
  }


  float pulseN[HcalConst::maxSamples];
  float pulseM[HcalConst::maxSamples];
  float pulseP[HcalConst::maxSamples];

  for (unsigned int i=0; i<HcalConst::maxSamples; i++) {
    pulseNn[idx+stride*i] = 0.f;
    pulseMn[idx+stride*i] = 0.f;
    pulsePn[idx+stride*i] = 0.f;
  }

 const float xx = t0;
 const float xxm = -nnlsWork_.dt + t0;
 const float xxp = nnlsWork_.dt + t0;

//  (*pfunctor_)(&xx[0]);
  functor_.singlePulseShapeFunc(&xx);
  functor_.getPulseShape(pulseNn, idx, stride);
  //  functor_.getPulseShape(pulseN);

//  (*pfunctor_)(&xxm[0]);
  functor_.singlePulseShapeFunc(&xxm);
  functor_.getPulseShape(pulseMn, idx, stride);
  //  functor_.getPulseShape(pulseM);
  
//  (*pfunctor_)(&xxp[0]);
  functor_.singlePulseShapeFunc(&xxp);
  functor_.getPulseShape(pulsePn, idx, stride);
  //  functor_.getPulseShape(pulseP);

  //in the 2018+ case where the sample of interest (SOI) is in TS3, add an extra offset to align 
  //with previous SOI=TS4 case assumed by psfPtr_->getPulseShape()
  int delta = 4 - nnlsWork_.tsOffset;

  auto invDt = 0.5f / nnlsWork_.dt;

  for (unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS) {

    pulseShape[iTS+nnlsWork_.maxoffset] = pulseNn[idx+stride*(iTS+delta)];
    if(calculateArrivalTime_) pulseDeriv[iTS+nnlsWork_.maxoffset] = (pulseMn[idx+stride*(iTS+delta)]-pulsePn[idx+stride*(iTS+delta)])*invDt;

    pulseMn[idx+stride*(iTS+delta)] -= pulseNn[idx+stride*(iTS+delta)];
    pulsePn[idx+stride*(iTS+delta)] -= pulseNn[idx+stride*(iTS+delta)];

  }

  for (unsigned int iTS=0; iTS<nnlsWork_.tsSize; ++iTS) {
    for (unsigned int jTS=0; jTS<iTS+1; ++jTS) {

      auto const  tmp = 0.5 * ( pulsePn[idx+stride*(iTS+delta)]*pulsePn[idx+stride*(jTS+delta)] + pulseMn[idx+stride*(iTS+delta)]*pulseMn[idx+stride*(jTS+delta)] );
      pulseCov(iTS+nnlsWork_.maxoffset,jTS+nnlsWork_.maxoffset) = tmp;

    }
  }
  
}

__device__
//void MahiFit::updateCov(const SampleMatrix& samplecov, SampleMatrixMAP const * pulseCovArray_) const {
void MahiFit::updateCov(SampleMatrixMAP const * pulseCovArray_) const {

  //  SampleMatrix invCovMat=samplecov;
  SampleMatrix invCovMat;
  invCovMat.setConstant(nnlsWork_.tsSize, nnlsWork_.tsSize, nnlsWork_.pedVal);
  invCovMat += nnlsWork_.noiseTerms.asDiagonal();

  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    auto const amp = nnlsWork_.ampVec.coeff(iBX);
    if (amp == 0) continue;

    int offset=nnlsWork_.bxs.coeff(iBX);

    if (offset==pedestalBX_) continue;		       
    else { 
      auto const ampsq = amp * amp;
      invCovMat += ampsq * pulseCovArray_[offset + nnlsWork_.bxOffset];
    }
  }

  nnlsWork_.covDecomp.compute(invCovMat);
}

__device__
float MahiFit::calculateArrivalTime(SamplePulseMatrixMAP & pulseMat_, SamplePulseMatrixMAP & pulseDerivMat_, unsigned int itIndex) const {

  if (nnlsWork_.nPulseTot > 1) {
    SamplePulseMatrix pulseDerivMatTMP = pulseDerivMat_;
    for (unsigned int iBX = 0; iBX < nnlsWork_.nPulseTot; ++iBX) {
      pulseDerivMat_.col(iBX) = pulseDerivMatTMP.col(nnlsWork_.bxs.coeff(iBX) + nnlsWork_.bxOffset);
    }
  }

  for (unsigned int iBX=0; iBX<nnlsWork_.nPulseTot; ++iBX) {
    pulseDerivMat_.col(iBX) *= nnlsWork_.ampVec.coeff(iBX);
  }

  SampleVector residuals = pulseMat_*nnlsWork_.ampVec - nnlsWork_.amplitudes;
  PulseVector solution = pulseDerivMat_.colPivHouseholderQr().solve(residuals);
  float t = solution.coeff(itIndex);
  t = (t>timeLimit_) ?  timeLimit_ : 
    ((t<-timeLimit_) ? -timeLimit_ : t);

  return t;

}
  
__device__
void MahiFit::nnls(SamplePulseMatrixMAP & pulseMat_) const {

  const unsigned int npulse = nnlsWork_.nPulseTot;
  const unsigned int nsamples = nnlsWork_.tsSize;

  nnlsWork_.invcovp = nnlsWork_.covDecomp.matrixL().solve(pulseMat_);
  nnlsWork_.aTaMat = nnlsWork_.invcovp.transpose().lazyProduct(nnlsWork_.invcovp);
  nnlsWork_.aTbVec = nnlsWork_.invcovp.transpose().lazyProduct(nnlsWork_.covDecomp.matrixL().solve(nnlsWork_.amplitudes));
  
  int iter = 0;
  Index idxwmax = 0;
  float wmax = 0.0;
  float threshold = nnlsThresh_;

  //  nnlsWork_.nP=0;
  
  while (true) {    
    if (iter>0 || nnlsWork_.nP==0) {
      if ( nnlsWork_.nP==std::min(npulse, nsamples)) break;
      
      const unsigned int nActive = npulse - nnlsWork_.nP;
      // exit if there are no more pulses to constrain
      if (nActive == 0)
        break;

      PulseVector updateWork = nnlsWork_.aTbVec - nnlsWork_.aTaMat*nnlsWork_.ampVec;
      
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
      nnlsUnconstrainParameter(idxp, pulseMat_);

    }

    while (true) {
      if (nnlsWork_.nP==0) break;     

      PulseVector ampvecpermtest = nnlsWork_.ampVec.head(nnlsWork_.nP);
      
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
      nnlsWork_.ampVec.coeffRef(minratioidx) = 0.f;
      
      nnlsConstrainParameter(minratioidx, pulseMat_);
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
void MahiFit::onePulseMinimize(SamplePulseMatrixMAP & pulseMat_) const {

  nnlsWork_.invcovp = nnlsWork_.covDecomp.matrixL().solve(pulseMat_);

  float aTaCoeff = (nnlsWork_.invcovp.transpose().lazyProduct(nnlsWork_.invcovp)).coeff(0);

  float aTbCoeff = nnlsWork_.invcovp.transpose().lazyProduct(nnlsWork_.covDecomp.matrixL().solve(nnlsWork_.amplitudes)).coeff(0);

  nnlsWork_.ampVec.coeffRef(0) = std::max(0.f, aTbCoeff / aTaCoeff);


}

__device__
float MahiFit::calculateChiSq(SamplePulseMatrixMAP & pulseMat_) const {
  
  return (nnlsWork_.covDecomp.matrixL().solve(pulseMat_*nnlsWork_.ampVec - nnlsWork_.amplitudes)).squaredNorm();
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
void MahiFit::nnlsUnconstrainParameter(Index idxp, SamplePulseMatrixMAP & pulseMat_) const {

  if (idxp != nnlsWork_.nP) {
    nnlsWork_.aTaMat.col(nnlsWork_.nP).swap(nnlsWork_.aTaMat.col(idxp));
    nnlsWork_.aTaMat.row(nnlsWork_.nP).swap(nnlsWork_.aTaMat.row(idxp));
    Eigen::numext::swap(nnlsWork_.aTbVec.coeffRef(nnlsWork_.nP),nnlsWork_.aTbVec.coeffRef(idxp));
    Eigen::numext::swap(nnlsWork_.ampVec.coeffRef(nnlsWork_.nP),nnlsWork_.ampVec.coeffRef(idxp)); // Victor did only this swap
    pulseMat_.col(nnlsWork_.nP).swap(pulseMat_.col(idxp));
    //  if(calculateArrivalTime_) pulseDerivMat_.col(nnlsWork_.nP).swap(pulseDerivMat_.col(idxp));
    Eigen::numext::swap(nnlsWork_.bxs.coeffRef(nnlsWork_.nP),nnlsWork_.bxs.coeffRef(idxp));
  }
  ++nnlsWork_.nP;

}

__device__
void MahiFit::nnlsConstrainParameter(Index minratioidx, SamplePulseMatrixMAP& pulseMat_) const {

  if (minratioidx != (nnlsWork_.nP - 1)) {
    nnlsWork_.aTaMat.col(nnlsWork_.nP-1).swap(nnlsWork_.aTaMat.col(minratioidx));
    nnlsWork_.aTaMat.row(nnlsWork_.nP-1).swap(nnlsWork_.aTaMat.row(minratioidx));
    Eigen::numext::swap(nnlsWork_.aTbVec.coeffRef(nnlsWork_.nP-1),nnlsWork_.aTbVec.coeffRef(minratioidx));
    Eigen::numext::swap(nnlsWork_.ampVec.coeffRef(nnlsWork_.nP-1),nnlsWork_.ampVec.coeffRef(minratioidx));
    pulseMat_.col(nnlsWork_.nP-1).swap(pulseMat_.col(minratioidx));
    //  if(calculateArrivalTime_) pulseDerivMat_.col(nnlsWork_.nP-1).swap(pulseDerivMat_.col(minratioidx));
    if (dynamicPed_ || calculateArrivalTime_) Eigen::numext::swap(nnlsWork_.bxs.coeffRef(nnlsWork_.nP-1),nnlsWork_.bxs.coeffRef(minratioidx));
  }
  --nnlsWork_.nP;

}

__device__
void MahiFit::solveSubmatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned nP) const {
  using namespace Eigen;
  switch( nP ) { // pulse matrix is always square.
    /*
  case 10:
    {
      Matrix<double,10,10> temp = mat;
      outvec.head<10>() = temp.ldlt().solve(invec.head<10>());
    }
    break;
  case 9:
    {
      Matrix<double,9,9> temp = mat.topLeftCorner<9,9>();
      outvec.head<9>() = temp.ldlt().solve(invec.head<9>());
    }
    break;
    */
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
//  nnlsWork_.fullTSOffset=0;
  nnlsWork_.bxOffset=0;
  nnlsWork_.maxoffset=0;
  nnlsWork_.dt=0;
  nnlsWork_.nP=0;

  nnlsWork_.amplitudes.setZero();
  nnlsWork_.noiseTerms.setZero();
//  nnlsWork_.pedConstraint.setZero();


}


}}
