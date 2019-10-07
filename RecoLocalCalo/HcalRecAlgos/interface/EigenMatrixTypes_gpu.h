#ifndef RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_h
#define RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_h

#include <Eigen/Dense>

constexpr int MaxSVSize=10;
constexpr int MaxFSVSize=15;
constexpr int MaxPVSize=8;

//constexpr int MaxSVSize=10;
//constexpr int MaxFSVSize=19;
//constexpr int MaxPVSize=10;

typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,MaxSVSize,1> SampleVector;
typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,MaxPVSize,1> PulseVector;
typedef Eigen::Matrix<int,Eigen::Dynamic,1,0,MaxPVSize,1> BXVector;

//typedef Eigen::Matrix<double,Eigen::Dynamic,1,0,MaxFSVSize,1> FullSampleVector;
//typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,MaxFSVSize,MaxFSVSize> FullSampleMatrix;
using FullSampleVector=Eigen::Matrix<double,MaxFSVSize,1>;
using FullSampleMatrix=Eigen::Matrix<double,MaxFSVSize,MaxFSVSize>;

typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,MaxSVSize,MaxSVSize> SampleMatrix;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,MaxPVSize,MaxPVSize> PulseMatrix;
typedef Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,0,MaxSVSize,MaxPVSize> SamplePulseMatrix;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, MaxPVSize, MaxSVSize> PulseSampleMatrix;

typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;
typedef Eigen::LLT<PulseMatrix> PulseDecompLLT;
typedef Eigen::LDLT<PulseMatrix> PulseDecompLDLT;

#endif
