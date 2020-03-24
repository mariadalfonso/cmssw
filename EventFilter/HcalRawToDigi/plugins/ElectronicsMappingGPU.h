#ifndef EventFilter_HcalRawToDigi_interface_ElectronicsMappingGPU_h
#define EventFilter_HcalRawToDigi_interface_ElectronicsMappingGPU_h

#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"

#ifndef __CUDACC__
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"
#include "HeterogeneousCore/CUDACore/interface/ESProduct.h"
#endif

namespace hcal { namespace raw {

class ElectronicsMappingGPU {
public:
    struct Product {
        ~Product();
        // trigger
        uint32_t *eid2tid;
        // detector
        uint32_t *eid2did;
    };

#ifndef __CUDACC__

    // rearrange pedestals
    ElectronicsMappingGPU(HcalElectronicsMap const&);

    // will call dealloation for Product thru ~Product
    ~ElectronicsMappingGPU() = default;

    // get device pointers
    Product const& getProduct(cudaStream_t) const;

    // 
    static std::string name() { return std::string{"hcalElectronicsMappingGPU"}; }

private:
    // in the future, we need to arrange so to avoid this copy on the host
    // if possible
    std::vector<uint32_t, CUDAHostAllocator<uint32_t>> eid2tid_;
    std::vector<uint32_t, CUDAHostAllocator<uint32_t>> eid2did_;

    cms::cuda::ESProduct<Product> product_;
#endif
};

}}

#endif // EventFilter_HcalRawToDigi_interface_ElectronicsMappingGPU_h
