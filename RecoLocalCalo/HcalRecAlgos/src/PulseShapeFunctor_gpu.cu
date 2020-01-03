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

    for (int i=0; i<HcalConst::maxPSshapeBin; i++) {
        acc25nsVec[i] = 0.f;
        diff25nsItvlVec[i] = 0.f;
    }

    for (int i=0; i<HcalConst::nsPerBX; i++) {
        accVarLenIdxZEROVec[i] = 0.f;
        diffVarItvlIdxZEROVec[i] = 0.f;
        accVarLenIdxMinusOneVec[i] = 0.f;
        diffVarItvlIdxMinusOneVec[i] = 0.f;
    }

    //The raw pulse
    for(int i=0;i<HcalConst::maxPSshapeBin;++i)  {
      pulse_hist[i] = pulse[i];
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
      // needed for M2
      //      psFit_x[i]      = 0;
      //      psFit_y[i]      = 0;
      //      psFit_erry[i]   = 1.;
      //      psFit_erry2[i]  = 1.;
      psFit_slew [i]  = 0.f;
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
  void PulseShapeFunctor::funcShape(float ntmpbin[HcalConst::maxSamples],
    const float pulseTime, /*const double pulseHeight,*/const float slew) {

    // pulse shape components over a range of time 0 ns to 255 ns in 1 ns steps
    constexpr int ns_per_bx = HcalConst::nsPerBX;
    //Get the starting time
    int i_start         = ( -HcalConst::iniTimeShift - pulseTime - slew >0 ? 0 : (int)std::abs(-HcalConst::iniTimeShift-pulseTime-slew) + 1);
    float offset_start = i_start - HcalConst::iniTimeShift - pulseTime - slew; //-199-2*pars[0]-2.*slew (for pars[0] > 98.5) or just -98.5-pars[0]-slew;
    // zeroing output binned pulse shape
    // TODO: do we need double here or float will be sufficient???
    for (unsigned int i=0; i<HcalConst::maxSamples; i++)
        ntmpbin[i] = 0.0f;

    if( edm::isNotFinite(offset_start) ){ //Check for nan
      ++ cntNANinfit;
    }else{
      if( offset_start == 1.0f ){ offset_start = 0.f; i_start-=1; } //Deal with boundary

      const int bin_start        = (int) offset_start; //bin off to integer
      const int bin_0_start      = ( offset_start < bin_start + 0.5f ? bin_start -1 : bin_start ); //Round it
      const int iTS_start        = i_start/ns_per_bx;         //Time Slice for time shift
      const int distTo25ns_start = ns_per_bx - 1 - i_start % ns_per_bx;    //Delta ns
      const float factor = offset_start - bin_0_start - 0.5f; //Small correction?
    
      //Build the new pulse
      ntmpbin[iTS_start] = (bin_0_start == -1 ? // Initial bin (I'm assuming this is ok)
			      accVarLenIdxMinusOneVec[distTo25ns_start] + factor * diffVarItvlIdxMinusOneVec[distTo25ns_start]
			    : accVarLenIdxZEROVec    [distTo25ns_start] + factor * diffVarItvlIdxZEROVec    [distTo25ns_start]);
      //Fill the rest of the bins
      for(int iTS = iTS_start+1; iTS < HcalConst::maxSamples; ++iTS){
	int bin_idx = distTo25ns_start + 1 + (iTS-iTS_start-1)*ns_per_bx + bin_0_start;
	ntmpbin[iTS] = acc25nsVec[bin_idx] + factor * diff25nsItvlVec[bin_idx];
      }
/*
      //Scale the pulse 
      for(int i=iTS_start; i < HcalConst::maxSamples; ++i) {
	ntmpbin[i]     *= pulseHeight;
      }
*/
    }

    return;
  }


  __device__
  void PulseShapeFunctor::EvalPulse(const float *pars) {

    int time = (pars[0]+timeShift_-timeMean_)*HcalConst::invertnsPerBx;
    funcShape(pulse_shape_, pars[0],psFit_slew[time]);

    return;

  }
  
}

}}
