#ifndef RecoLocalCalo_HcalRecAlgos_interface_HcalQIECodersGPU_h
#define RecoLocalCalo_HcalRecAlgos_interface_HcalQIECodersGPU_h

#include "CondFormats/HcalObjects/interface/HcalQIEData.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/CUDAESProduct.h"
#endif

class HcalQIECodersGPU {
public:
    static constexpr uint32_t numValuesPerChannel = 16;

    struct Product {
        ~Product();
        float *offsets;
        float *slopes;
    };

#ifndef __CUDACC__
    // rearrange reco params
    HcalQIECodersGPU(HcalQIEData const&);

    // will trigger deallocation of Product thru ~Product
    ~HcalQIECodersGPU() = default;

    // get device pointers
    Product const& getProduct(cuda::stream_t<>&) const;

    // 
    static std::string name() { return std::string{"hcalQIECodersGPU"}; }

private:
    uint64_t totalChannels_;
    std::vector<float, CUDAHostAllocator<float>> offsets_;
    std::vector<float, CUDAHostAllocator<float>> slopes_;

    CUDAESProduct<Product> product_;
#endif
};

#endif
