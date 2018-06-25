#ifndef RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_mahi_h
#define RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_mahi_h

#include <vector>

#include <cuda.h>
#include <cuda_runtime.h>

#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalCalibrations.h"
#include "CondFormats/HcalObjects/interface/HcalRecoParam.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/EigenMatrixTypes.h"

namespace hcal { namespace mahi {

struct Workspace {
    unsigned int nPulseTot;
    unsigned int tsSize;
    unsigned int tsOffset;
    unsigned int fullTSOffset;
    int bxOffset;
    double dt;

    BXVector bxs;
    SampleVector amplitudes;
    SampleMatrix invCovMat;
    SampleVector noiseTerms;
    SampleMatrix pedConstraint;

    FullSampleMatrix pulseCovArray[MaxPVSize];
    FullSampleVector pulseShapeArray[MaxPVSize];
    FullSampleVector pulseDerivArray[MaxPVSize];

    double pulseN[MaxSVSize];
    double pulseM[MaxSVSize];
    double pulseP[MaxSVSize];

    SamplePulseMatrix pulseMat;
    SamplePulseMatrix pulseDerivMat;
    PulseVector residuals;

    unsigned int nP;
    PulseVector ampVec;
    PulseVector errVec;
    PulseVector ampvecpermtest;

    SamplePulseMatrix invcovp;
    PulseMatrix aTaMat;
    PulseVector aTbVec;
    PulseVector updateWork;

    // TODO: solvers, etc...
};

struct RecValues {
    float energy;
    float time;
    float chi2;
};

struct DeviceData {
    HBHEChannelInfo         *vinfos;
    HBHERecHit              *vrechits;
    HcalRecoParam           *vparams;
    HcalCalibrations        *vcalibs;

    void allocate(int size) {
        cudaMalloc((void**)&vinfos, size * sizeof(HBHEChannelInfo));
        cudaMalloc((void**)&vrechits, size * sizeof(HBHERecHit));
        cudaMalloc((void**)&vparams, size * sizeof(HcalRecoParam));
        cudaMalloc((void**)&vcalibs, size* sizeof(HcalCalibrations));
    }
    void free() {
        cudaFree(vinfos);
        cudaFree(vrechits);
        cudaFree(vparams);
        cudaFree(vcalibs);
    }
};

// allocate on the de

// reconstruction
void reco(DeviceData,
          HBHEChannelInfoCollection&, HBHERecHitCollection&, 
          std::vector<HcalRecoParam> const&, std::vector<HcalCalibrations> const&, bool);

}}

#endif // RecoLocalCalo_HcalRecAlgos_interface_gpu_reco_mahi_h
