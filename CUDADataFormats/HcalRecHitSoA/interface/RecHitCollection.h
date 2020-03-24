#ifndef CUDADataFormats_HcalRecHitCollectionSoA_interface_RecHitCollection_h
#define CUDADataFormats_HcalRecHitCollectionSoA_interface_RecHitCollection_h

#include <vector>

#include "CUDADataFormats/HcalCommon/interface/Common.h"
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"

namespace hcal {

template<typename StoragePolicy>
struct RecHitCollection 
        : public common::AddSize<typename StoragePolicy::TagType> {
    RecHitCollection() = default;
    RecHitCollection(const RecHitCollection&) = default;
    RecHitCollection& operator=(const RecHitCollection&) = default;

    RecHitCollection(RecHitCollection&&) = default;
    RecHitCollection& operator=(RecHitCollection&&) = default;

    typename StoragePolicy::template StorageSelector<float>::type energy;
    typename StoragePolicy::template StorageSelector<float>::type chi2;
    typename StoragePolicy::template StorageSelector<float>::type energyM0;
    typename StoragePolicy::template StorageSelector<float>::type timeM0;
    typename StoragePolicy::template StorageSelector<uint32_t>::type did;

    template<typename U = typename StoragePolicy::TagType>
    typename std::enable_if<std::is_same<U, common::tags::Vec>::value, void>::type 
    resize(size_t size) {
        energy.resize(size);
        chi2.resize(size);
        energyM0.resize(size);
        timeM0.resize(size);
        did.resize(size);
    }
};

}

#endif // RecoLocalCalo_HcalRecAlgos_interface_RecHitCollection_h
