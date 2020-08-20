#ifndef DataFormats_L1TCalorimeter_HGCalTower_h
#define DataFormats_L1TCalorimeter_HGCalTower_h

#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1THGCal/interface/HGCalTowerID.h"
#include "DataFormats/L1THGCal/interface/HFNoseTowerID.h"

namespace l1t {

  template <typename T>
  class HGCalTowerT : public L1Candidate {
  public:
    HGCalTowerT() : etEm_(0.), etHad_(0.), id_(0), hwEtEm_(0), hwEtHad_(0), hwEtRatio_(0) {}

    HGCalTowerT(double etEm,
               double etHad,
               double eta,
               double phi,
               unsigned short id,
               int hwpt = 0,
               int hweta = 0,
               int hwphi = 0,
               int qual = 0,
               int hwEtEm = 0,
               int hwEtHad = 0,
               int hwEtRatio = 0);

    ~HGCalTowerT() override;

    void addEtEm(double et);
    void addEtHad(double et);

    double etEm() const { return etEm_; };
    double etHad() const { return etHad_; };

    const HGCalTowerT<T>& operator+=(const HGCalTowerT<T> & tower);

    T id() const { return id_; }
    short zside() const { return id_.zside(); }

    void setHwEtEm(int et) { hwEtEm_ = et; }
    void setHwEtHad(int et) { hwEtHad_ = et; }
    void setHwEtRatio(int ratio) { hwEtRatio_ = ratio; }

    int hwEtEm() const { return hwEtEm_; }
    int hwEtHad() const { return hwEtHad_; }
    int hwEtRatio() const { return hwEtRatio_; }

  private:
    void addEt(double et);

    // additional hardware quantities
    double etEm_;
    double etHad_;
    T id_;

    int hwEtEm_;
    int hwEtHad_;
    int hwEtRatio_;
  };

  using HGCalTower = HGCalTowerT<HGCalTowerID>;
  using HFNoseTower = HGCalTowerT<HFNoseTowerID>;

  typedef BXVector<HGCalTower> HGCalTowerBxCollection;
  typedef BXVector<HFNoseTower> HFNoseTowerBxCollection;

}  // namespace l1t

#endif
