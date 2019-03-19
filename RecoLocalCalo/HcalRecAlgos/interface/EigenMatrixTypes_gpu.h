#ifndef RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_gpu_h
#define RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_gpu_h

#include <Eigen/Dense>

constexpr int MaxSVSize=8;
constexpr int MaxFSVSize=19;
constexpr int MaxPVSize=10;

using SampleVector=Eigen::Matrix<double,MaxSVSize,1>;
typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,MaxPVSize,1> PulseVector;
//using PulseVector=Eigen::Matrix<double,MaxPVSize,1>;

typedef Eigen::Matrix<int,Eigen::Dynamic,1,0,MaxPVSize,1> BXVector;
//using BXVector=Eigen::Matrix<double,MaxPVSize,1>;

using FullSampleVector=Eigen::Matrix<double,MaxFSVSize,1>;
using FullSampleMatrix=Eigen::Matrix<double,MaxFSVSize,MaxFSVSize>;

using SampleMatrix=Eigen::Matrix<double,MaxSVSize,MaxSVSize>;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,MaxPVSize,MaxPVSize> PulseMatrix;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,MaxSVSize,MaxPVSize> SamplePulseMatrix;
//using PulseMatrix=Eigen::Matrix<double,MaxPVSize,MaxPVSize>;
//using SamplePulseMatrix=Eigen::Matrix<double,MaxSVSize,MaxPVSize>;

typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;
typedef Eigen::LLT<PulseMatrix> PulseDecompLLT;
typedef Eigen::LDLT<PulseMatrix> PulseDecompLDLT;

#endif
