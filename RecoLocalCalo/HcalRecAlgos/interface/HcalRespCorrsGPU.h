#ifndef RecoLocalCalo_HcalRecAlgos_interface_HcalRespCorrsGPU_h
#define RecoLocalCalo_HcalRecAlgos_interface_HcalRespCorrsGPU_h

#include "CondFormats/HcalObjects/interface/HcalRespCorrs.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/ESProduct.h"
#endif

class HcalRespCorrsGPU {
public:
    struct Product {
        ~Product();
        float *values;
    };

#ifndef __CUDACC__
    // rearrange reco params
    HcalRespCorrsGPU(HcalRespCorrs const&);

    // will trigger deallocation of Product thru ~Product
    ~HcalRespCorrsGPU() = default;

    // get device pointers
    Product const& getProduct(cudaStream_t) const;

    // 
    static std::string name() { return std::string{"hcalRespCorrsGPU"}; }

private:
    std::vector<float, CUDAHostAllocator<float>> values_;

    cms::cuda::ESProduct<Product> product_;
#endif
};

#endif
