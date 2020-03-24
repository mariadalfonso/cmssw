#ifndef RecoLocalCalo_HcalRecAlgos_interface_HcalSiPMCharacteristicsGPU_h
#define RecoLocalCalo_HcalRecAlgos_interface_HcalSiPMCharacteristicsGPU_h

#include "CondFormats/HcalObjects/interface/HcalSiPMCharacteristics.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

class HcalSiPMCharacteristicsGPU {
public:
    struct Product {
        ~Product();
        int *pixels;
        float *parLin1, *parLin2, *parLin3;
        float *crossTalk;
        int *auxi1;
        float *auxi2;
    };

#ifndef __CUDACC__
    // rearrange reco params
    HcalSiPMCharacteristicsGPU(HcalSiPMCharacteristics const&);

    // will trigger deallocation of Product thru ~Product
    ~HcalSiPMCharacteristicsGPU() = default;

    // get device pointers
    Product const& getProduct(cuda::stream_t<>&) const;

    // 
    static std::string name() { return std::string{"hcalSiPMCharacteristicsGPU"}; }

private:
    std::vector<int, CUDAHostAllocator<int>> pixels_, auxi1_;
    std::vector<float, CUDAHostAllocator<float>> parLin1_, parLin2_, parLin3_,
        crossTalk_, auxi2_;

    CUDAESProduct<Product> product_;
#endif
};

#endif
