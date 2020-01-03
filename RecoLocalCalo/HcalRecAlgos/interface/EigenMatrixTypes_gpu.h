#ifndef RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_h
#define RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_h

#include <Eigen/Dense>

constexpr int MaxSVSize=8;
constexpr int MaxFSVSize=15;
constexpr int MaxPVSize=8;

typedef Eigen::Matrix<float,Eigen::Dynamic,1,0,MaxSVSize,1> SampleVector;
typedef Eigen::Matrix<float,Eigen::Dynamic,1,0,MaxPVSize,1> PulseVector;
typedef Eigen::Matrix<int,Eigen::Dynamic,1,0,MaxPVSize,1> BXVector;

//typedef Eigen::Matrix<float,Eigen::Dynamic,1,0,MaxFSVSize,1> FullSampleVector;
//typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,0,MaxFSVSize,MaxFSVSize> FullSampleMatrix;
using FullSampleVector=Eigen::Matrix<float,MaxFSVSize,1>;
using FullSampleMatrix=Eigen::Matrix<float,MaxFSVSize,MaxFSVSize>;

typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,0,MaxSVSize,MaxSVSize> SampleMatrix;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,0,MaxPVSize,MaxPVSize> PulseMatrix;
typedef Eigen::Matrix<float,Eigen::Dynamic,Eigen::Dynamic,0,MaxSVSize,MaxPVSize> SamplePulseMatrix;

typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, 0, MaxPVSize, MaxSVSize> PulseSampleMatrix;

typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;

/////

using SampleMatrixFixed=Eigen::Matrix<float,MaxSVSize,MaxSVSize>;
using SamplePulseMatrixFixed =Eigen::Matrix<float,MaxSVSize,MaxPVSize>;

using DynStride = Eigen::Stride<Eigen::Dynamic,Eigen::Dynamic>;
using SamplePulseMatrixMAP = Eigen::Map<SamplePulseMatrixFixed, 0, DynStride>;
using SampleMatrixMAP = Eigen::Map<SampleMatrixFixed, 0, DynStride>;

#endif
