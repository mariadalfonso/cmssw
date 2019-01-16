#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFunctor_gpu.h"
#include "FWCore/Utilities/interface/isFinite.h"

namespace hcal { namespace mahi {

namespace FitterFuncs{

  //Decalare the Pulse object take it in from Hcal and set some options
  __device__
  PulseShapeFunctor::PulseShapeFunctor(float const* pulse,
				                       bool iPedestalConstraint, 
                                       bool iTimeConstraint,bool iAddPulseJitter,
				                       double iPulseJitter,double iTimeMean,
                                       double iPedMean,
				                       unsigned nSamplesToFit) {
    cntNANinfit = 0;

    //The raw pulse
    for(int i=0;i<HcalConst::maxPSshapeBin;++i)  {
        // done according to HcalPulseShape.cc::at(double t)
        pulse_hist[i] = pulse[(int)(static_cast<double>(i) + 0.5)];
    }

    // Accumulate 25ns for each starting point of 0, 1, 2, 3...
    for(int i=0; i<HcalConst::maxPSshapeBin; ++i){
      for(int j=i; j<i+HcalConst::nsPerBX; ++j){  //sum over HcalConst::nsPerBXns from point i
	acc25nsVec[i] += ( j < HcalConst::maxPSshapeBin? pulse_hist[j] : pulse_hist[HcalConst::maxPSshapeBin-1]);
      }
      diff25nsItvlVec[i] = ( i+HcalConst::nsPerBX < HcalConst::maxPSshapeBin? pulse_hist[i+HcalConst::nsPerBX] - pulse_hist[i] : pulse_hist[HcalConst::maxPSshapeBin-1] - pulse_hist[i]);
    }
    // Accumulate different ns for starting point of index either 0 or -1
    for(int i=0; i<HcalConst::nsPerBX; ++i){
      if( i==0 ){
	accVarLenIdxZEROVec[0] = pulse_hist[0];
	accVarLenIdxMinusOneVec[i] = pulse_hist[0];
      } else{
	accVarLenIdxZEROVec[i] = accVarLenIdxZEROVec[i-1] + pulse_hist[i];
	accVarLenIdxMinusOneVec[i] = accVarLenIdxMinusOneVec[i-1] + pulse_hist[i-1];
      }
      diffVarItvlIdxZEROVec[i] = pulse_hist[i+1] - pulse_hist[0];
      diffVarItvlIdxMinusOneVec[i] = pulse_hist[i] - pulse_hist[0];
    }
    for(int i = 0; i < HcalConst::maxSamples; i++) { 
      psFit_x[i]      = 0;
      psFit_y[i]      = 0;
      psFit_erry[i]   = 1.;
      psFit_erry2[i]  = 1.;
      psFit_slew [i]  = 0.;
    }
    //Constraints
    pedestalConstraint_ = iPedestalConstraint;
    timeConstraint_     = iTimeConstraint;
    addPulseJitter_     = iAddPulseJitter;
    pulseJitter_        = iPulseJitter*iPulseJitter;

    // for M2
    timeMean_           = iTimeMean;
    pedMean_            = iPedMean;
    timeShift_          = 100.;
    timeShift_ += 12.5;//we are trying to get BX 

    nSamplesToFit_ = nSamplesToFit;

  }

  __device__
  void PulseShapeFunctor::funcShape(double ntmpbin[HcalConst::maxSamples],
    const double pulseTime, const double pulseHeight,const double slew) {

    // pulse shape components over a range of time 0 ns to 255 ns in 1 ns steps
    constexpr int ns_per_bx = HcalConst::nsPerBX;
    constexpr int num_ns = HcalConst::nsPerBX*HcalConst::maxSamples;
    constexpr int num_bx = num_ns/ns_per_bx;
    //Get the starting time
    int i_start         = ( -HcalConst::iniTimeShift - pulseTime - slew >0 ? 0 : (int)std::abs(-HcalConst::iniTimeShift-pulseTime-slew) + 1);
    double offset_start = i_start - HcalConst::iniTimeShift - pulseTime - slew; //-199-2*pars[0]-2.*slew (for pars[0] > 98.5) or just -98.5-pars[0]-slew;
    // zeroing output binned pulse shape
    // TODO: do we need double here or float will be sufficient???
    for (unsigned int i=0; i<HcalConst::maxSamples; i++)
        ntmpbin[i] = 0.0;

    if( edm::isNotFinite(offset_start) ){ //Check for nan
      ++ cntNANinfit;
    }else{
      if( offset_start == 1.0 ){ offset_start = 0.; i_start-=1; } //Deal with boundary

      const int bin_start        = (int) offset_start; //bin off to integer
      const int bin_0_start      = ( offset_start < bin_start + 0.5 ? bin_start -1 : bin_start ); //Round it
      const int iTS_start        = i_start/ns_per_bx;         //Time Slice for time shift
      const int distTo25ns_start = HcalConst::nsPerBX - 1 - i_start%ns_per_bx;    //Delta ns 
      const double factor = offset_start - bin_0_start - 0.5; //Small correction?
    
      //Build the new pulse
      ntmpbin[iTS_start] = (bin_0_start == -1 ? // Initial bin (I'm assuming this is ok)
			      accVarLenIdxMinusOneVec[distTo25ns_start] + factor * diffVarItvlIdxMinusOneVec[distTo25ns_start]
			    : accVarLenIdxZEROVec    [distTo25ns_start] + factor * diffVarItvlIdxZEROVec    [distTo25ns_start]);
      //Fill the rest of the bins
      for(int iTS = iTS_start+1; iTS < num_bx; ++iTS){
	int bin_idx = distTo25ns_start + 1 + (iTS-iTS_start-1)*ns_per_bx + bin_0_start;
	ntmpbin[iTS] = acc25nsVec[bin_idx] + factor * diff25nsItvlVec[bin_idx];
      }
      //Scale the pulse 
      for(int i=iTS_start; i < num_bx; ++i) {
	ntmpbin[i]     *= pulseHeight;
      }

    }

    return;
  }

  __device__
  double PulseShapeFunctor::EvalPulse(const double *pars, const unsigned nPars) {

      unsigned i =0, j=0;

      const double pedestal=pars[nPars-1];

      //Stop crashes
      for(i =0; i < nPars; ++i ) if( edm::isNotFinite(pars[i]) ){ ++ cntNANinfit; return 1e10; }
      
      //calculate chisquare
      double chisq  = 0;
      const unsigned parBy2=(nPars-1)/2;
      //      std::array<float,HcalConst::maxSamples> pulse_shape_;

      if(addPulseJitter_) {
	int time = (pars[0]+timeShift_-timeMean_)*HcalConst::invertnsPerBx;
	//Interpolate the fit (Quickly)
	funcShape(pulse_shape_, pars[0],pars[1],psFit_slew[time]);
	for (j=0; j<nSamplesToFit_; ++j) {
	  psFit_erry2[j]  += pulse_shape_[j]*pulse_shape_[j]*pulseJitter_;
	  pulse_shape_sum_[j] = pulse_shape_[j] + pedestal;
	}

	for (i=1; i<parBy2; ++i) {
	  time = (pars[i*2]+timeShift_-timeMean_)*HcalConst::invertnsPerBx;
	  //Interpolate the fit (Quickly)
	  funcShape(pulse_shape_, pars[i*2],pars[i*2+1],psFit_slew[time]);
	  // add an uncertainty from the pulse (currently noise * pulse height =>Ecal uses full cov)
	 /////
	  for (j=0; j<nSamplesToFit_; ++j) {
	    psFit_erry2[j] += pulse_shape_[j]*pulse_shape_[j]*pulseJitter_;
	    pulse_shape_sum_[j] += pulse_shape_[j];
	  }	    
	}
      }
      else{
	int time = (pars[0]+timeShift_-timeMean_)*HcalConst::invertnsPerBx;
	//Interpolate the fit (Quickly)
	funcShape(pulse_shape_, pars[0],pars[1],psFit_slew[time]);
	for(j=0; j<nSamplesToFit_; ++j)
	  pulse_shape_sum_[j] = pulse_shape_[j] + pedestal;

	for (i=1; i<parBy2; ++i) {
	  time = (pars[i*2]+timeShift_-timeMean_)*HcalConst::invertnsPerBx;
	  //Interpolate the fit (Quickly)
	  funcShape(pulse_shape_, pars[i*2],pars[i*2+1],psFit_slew[time]);
	  // add an uncertainty from the pulse (currently noise * pulse height =>Ecal uses full cov)
	  for(j=0; j<nSamplesToFit_; ++j)
	    pulse_shape_sum_[j] += pulse_shape_[j];
	}
      }

      for (i=0;i<nSamplesToFit_; ++i)
	{
	  const double d = psFit_y[i]- pulse_shape_sum_[i];
	  chisq += d*d/psFit_erry2[i];
	}

      if(pedestalConstraint_) {
	 //Add the pedestal Constraint to chi2
         chisq += invertpedSig2_*(pedestal - pedMean_)*(pedestal - pedMean_);
      }
        //Add the time Constraint to chi2
      if(timeConstraint_) {
	for(j=0; j< parBy2; ++j ){
	  int time = (pars[j*2]+timeShift_-timeMean_)*(double)HcalConst::invertnsPerBx;
	  double time1 = -100.+time*HcalConst::nsPerBX;
	  chisq += inverttimeSig2_*(pars[j*2] - timeMean_ - time1)*(pars[j*2] - timeMean_ - time1);
	}
      }
      return chisq;
   }

  __device__
  double PulseShapeFunctor::singlePulseShapeFunc( const double *x ) {
    return EvalPulse(x,3);
  }
  
  __device__
  double PulseShapeFunctor::doublePulseShapeFunc( const double *x ) {
    return EvalPulse(x,5);
  }
  
  __device__
  double PulseShapeFunctor::triplePulseShapeFunc( const double *x ) {
    return EvalPulse(x,7);
  }

}

}}
