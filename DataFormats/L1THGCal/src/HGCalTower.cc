#include "DataFormats/L1THGCal/interface/HGCalTower.h"
#include "FWCore/Utilities/interface/EDMException.h"

template class l1t::HGCalTowerT<l1t::HGCalTowerID>;
template class l1t::HGCalTowerT<l1t::HFNoseTowerID>;
using l1t::L1Candidate;

template <typename T>
l1t::HGCalTowerT<T>::HGCalTowerT(double etEm,
                       double etHad,
                       double eta,
                       double phi,
                       unsigned short id,
                       int hwpt,
                       int hweta,
                       int hwphi,
                       int qual,
                       int hwEtEm,
                       int hwEtHad,
                       int hwEtRatio)
    : L1Candidate(PolarLorentzVector(etEm + etHad, eta, phi, 0.), hwpt, hweta, hwphi, qual),
      etEm_(etEm),
      etHad_(etHad),
      id_(id),
      hwEtEm_(hwEtEm),
      hwEtHad_(hwEtHad),
      hwEtRatio_(hwEtRatio) {}

template <typename T>
l1t::HGCalTowerT<T>::~HGCalTowerT() {}

template <typename T>
void l1t::HGCalTowerT<T>::addEtEm(double et) {
  etEm_ += et;
  addEt(et);
}

template <typename T>
void l1t::HGCalTowerT<T>::addEtHad(double et) {
  etHad_ += et;
  addEt(et);
}

template <typename T>
void l1t::HGCalTowerT<T>::addEt(double et) { this->setP4(PolarLorentzVector(this->pt() + et, this->eta(), this->phi(), 0.)); }

template <typename T>
const l1t::HGCalTowerT<T>& l1t::HGCalTowerT<T>::operator+=(const HGCalTowerT<T>& tower) {
  // NOTE: assume same eta and phi -> need an explicit check on the ID
  if (id().rawId() != tower.id().rawId()) {
    throw edm::Exception(edm::errors::StdException, "StdException")
        << "HGCalTower: adding to this tower with ID: " << id().rawId()
        << " one with different ID: " << tower.id().rawId() << std::endl;
  }
  addEt(tower.pt());
  etEm_ += tower.etEm();
  etHad_ += tower.etHad();

  return *this;
}

