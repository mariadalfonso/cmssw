#ifndef RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_h
#define RecoLocalCalo_HcalRecAlgos_EigenMatrixTypes_h

#include <Eigen/Dense>

constexpr int MaxSVSize = 8;
constexpr int MaxFSVSize = 15;
constexpr int MaxPVSize = 8;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MaxSVSize, 1> SampleVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MaxPVSize, 1> PulseVector;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1, 0, MaxPVSize, 1> BXVector;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1, 0, MaxFSVSize, 1> FullSampleVector;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, MaxFSVSize, MaxFSVSize> FullSampleMatrix;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, MaxSVSize, MaxSVSize> SampleMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, MaxPVSize, MaxPVSize> PulseMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, MaxSVSize, MaxPVSize> SamplePulseMatrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, 0, MaxPVSize, MaxSVSize> PulseSampleMatrix;

typedef Eigen::LLT<SampleMatrix> SampleDecompLLT;
typedef Eigen::LLT<PulseMatrix> PulseDecompLLT;

#endif
