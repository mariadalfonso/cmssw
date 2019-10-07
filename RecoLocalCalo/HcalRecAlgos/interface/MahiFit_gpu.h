#ifndef RecoLocalCalo_HcalRecAlgos_MahiFit_gpu_HH
#define RecoLocalCalo_HcalRecAlgos_MahiFit_gpu_HH

#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes_gpu.h"
#include "DataFormats/HcalRecHit/interface/HBHEChannelInfo.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFunctor_gpu.h"

namespace hcal { namespace mahi {

struct MahiNnlsWorkspace {

  unsigned int nPulseTot;
  unsigned int tsSize;
  unsigned int tsOffset;
  int bxOffset;
  int maxoffset;
  float dt;

  //holds active bunch crossings
  BXVector bxs;  

  //holds data samples
  SampleVector amplitudes;

  //holds diagonal noise terms
  SampleVector noiseTerms;

  //holds flat pedestal uncertainty
  float pedVal;
  
  //holds full covariance matrix for a pulse shape 
  //varied in time
  //  SampleMatrix pulseCovArray[MaxPVSize];

  //holds matrix of pulse shape templates for each BX
  //  SamplePulseMatrix pulseMat;

  //holds matrix of pulse shape derivatives for each BX
  //  SamplePulseMatrix pulseDerivMat;

  //for FNNLS algorithm
  unsigned int nP;
  PulseVector ampVec;

  SamplePulseMatrix invcovp;
  PulseMatrix aTaMat; // A-transpose A (matrix)
  PulseVector aTbVec; // A-transpose b (vector)

  SampleDecompLLT covDecomp;
  PulseDecompLDLT pulseDecomp;

};

struct MahiDebugInfo {

  int   nSamples;
  int   soi;

  bool  use3;

  float inTimeConst;
  float inDarkCurrent;
  float inPedAvg;
  float inGain;
  
  float inNoiseADC[MaxSVSize];
  float inNoiseDC[MaxSVSize];
  float inNoisePhoto[MaxSVSize];
  float inPedestal[MaxSVSize];

  float totalUCNoise[MaxSVSize];

  float mahiEnergy;
  float chiSq;
  float arrivalTime;

  float pEnergy;
  float nEnergy;
  float pedEnergy;

  float count[MaxSVSize];
  float inputTS[MaxSVSize];
  int inputTDC[MaxSVSize];
  float itPulse[MaxSVSize];
  float pPulse[MaxSVSize];
  float nPulse[MaxSVSize];
  
};

class MahiFit
{
 public:
  __device__ 
  MahiFit(float const*);

  __device__
  void phase1Apply(const HBHEChannelInfo& channelData, 
		   float& reconstructedEnergy, 
		   float& reconstructedTime, 
		   bool& useTriple,
		   float& chi2,
		   float*, float*, float*,
		   double*, double*, double*
		   ) const;

  __device__
  void phase1Debug(const HBHEChannelInfo& channelData,
		   MahiDebugInfo& mdi) const;

  __device__
    void doFit(float correctedOutput[3], const int nbx, float*, float*, float*,double*, double*, double*) const;

  __device__
  //void setPulseShapeTemplate  (const HcalPulseShapes::Shape& ps,const HcalTimeSlew * hcalTimeSlewDelay);
  void setPulseShapeTemplate(float const* pshape);
/*  __device__
  void resetPulseShapeTemplate(const HcalPulseShapes::Shape& ps);
  */

  typedef BXVector::Index Index;
  float const* pshape_;
  // TODO: resolve this
//  const HcalTimeSlew* hcalTimeSlewDelay_=nullptr;

 private:

  __device__
  double minimize(SamplePulseMatrixMAP & pulseMat_, SamplePulseMatrixMAP & pulseDerivMat_, SampleMatrixMAP const * pulseCovArray) const;
  __device__
  void onePulseMinimize(SamplePulseMatrixMAP & pulseMat_) const;
  __device__
  void updateCov(SampleMatrixMAP const * pulseCovArray) const;
  __device__
  void updatePulseShape(double itQ, FullSampleVector &pulseShape, 
			FullSampleVector &pulseDeriv,
			FullSampleMatrix &pulseCov,
			float*, float*, float*
			) const;

  __device__
  float calculateArrivalTime(SamplePulseMatrixMAP & pulseMat_, SamplePulseMatrixMAP & pulseDerivMat_) const;
  __device__
  double calculateChiSq(SamplePulseMatrixMAP & pulseMat_) const;
  __device__
  void nnls(SamplePulseMatrixMAP & pulseMat_, SamplePulseMatrixMAP & pulseDerivMat_) const;
  __device__
  void resetWorkspace() const;

  __device__
  void nnlsUnconstrainParameter(Index idxp, SamplePulseMatrixMAP & pulseMat_, SamplePulseMatrixMAP & pulseDerivMat_) const;
  __device__
  void nnlsConstrainParameter(Index minratioidx, SamplePulseMatrixMAP & pulseMat_, SamplePulseMatrixMAP & pulseDerivMat_) const;

  __device__
  void solveSubmatrix(PulseMatrix& mat, PulseVector& invec, PulseVector& outvec, unsigned nP) const;

  mutable MahiNnlsWorkspace nnlsWork_;

  //hard coded in initializer
  static constexpr int pedestalBX_ = 100;

  // used to restrict returned time value to a 25 ns window centered 
  // on the nominal arrival time
  static constexpr float timeLimit_ = 12.5;

  bool calculateArrivalTime_{false};

  // Python-configurables
  bool dynamicPed_ {false};
  float ts4Thresh_ {0.0};
  //  float chiSqSwitch_{15.0};
  float chiSqSwitch_{-1.};

  bool applyTimeSlew_{false};
  /* HcalTimeSlew::BiasSetting */ int slewFlavor_{1};
  float tsDelay1GeV_{10.0};

  float meanTime_{0.0};
  float timeSigmaHPD_{5.0}; 
  float timeSigmaSiPM_{2.5};

  int activeBXs_[8] = {-3, -2, -1, 0, 1, 2, 3, 4};

  int nMaxItersMin_{500}; 
  int nMaxItersNNLS_{500}; 

  float deltaChiSqThresh_{1e-3}; 
  float nnlsThresh_{1e-11}; 

  unsigned int bxSizeConf_{8};
  int bxOffsetConf_{3}; // -(min_element({-3, -2, -1, 0, 1, 2, 3, 4}))

  //for pulse shapes
  int cntsetPulseShape_;

  mutable FitterFuncs::PulseShapeFunctor functor_;
};

}}

#endif
