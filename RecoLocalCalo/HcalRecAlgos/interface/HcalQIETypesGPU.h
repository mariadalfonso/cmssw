#ifndef RecoLocalCalo_HcalRecAlgos_interface_HcalQIETypesGPU_h
#define RecoLocalCalo_HcalRecAlgos_interface_HcalQIETypesGPU_h

#include "CondFormats/HcalObjects/interface/HcalQIETypes.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

class HcalQIETypesGPU {
public:
    struct Product {
        ~Product();
        int *values;
    };

#ifndef __CUDACC__
    // rearrange reco params
    HcalQIETypesGPU(HcalQIETypes const&);

    // will trigger deallocation of Product thru ~Product
    ~HcalQIETypesGPU() = default;

    // get device pointers
    Product const& getProduct(cudaStream_t) const;

    // 
    static std::string name() { return std::string{"hcalQIETypesGPU"}; }

private:
    std::vector<int, CUDAHostAllocator<int>> values_;

    CUDAESProduct<Product> product_;
#endif
};

#endif
