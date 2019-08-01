#include "DataFormats/ForwardDetId/interface/HFNoseTriggerDetId.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <ostream>
#include <iostream>

const HFNoseTriggerDetId HFNoseTriggerDetId::Undefined(HFNose,0,0,0,0,0,0,0);

HFNoseTriggerDetId::HFNoseTriggerDetId() : DetId() {
}

HFNoseTriggerDetId::HFNoseTriggerDetId(uint32_t rawid) : DetId(rawid) {
}

HFNoseTriggerDetId::HFNoseTriggerDetId(int subdet, int zp, int type, int layer,
				     int waferU, int waferV, int cellU,
				     int cellV) : DetId(HGCalTrigger,ForwardEmpty) {  

  int waferUabs(std::abs(waferU)), waferVabs(std::abs(waferV));
  int waferUsign = (waferU >= 0) ? 0 : 1;
  int waferVsign = (waferV >= 0) ? 0 : 1;
  int zside      = (zp < 0) ? 1 : 0;
  int lay   = std::max(layer-1, 0);
  id_ |= (((cellU     & kHFNoseCellUMask) << kHFNoseCellUOffset) |
	  ((cellV     & kHFNoseCellVMask) << kHFNoseCellVOffset) |
	  ((waferUabs & kHFNoseWaferUMask) << kHFNoseWaferUOffset) |
	  ((waferUsign& kHFNoseWaferUSignMask) << kHFNoseWaferUSignOffset) |
	  ((waferVabs & kHFNoseWaferVMask) << kHFNoseWaferVOffset) |
	  ((waferVsign& kHFNoseWaferVSignMask) << kHFNoseWaferVSignOffset) |
	  ((lay     & kHFNoseLayerMask) << kHFNoseLayerOffset) |
	  ((zside     & kHFNoseZsideMask) << kHFNoseZsideOffset) |
	  ((type      & kHFNoseTypeMask)  << kHFNoseTypeOffset) |
	  ((subdet    & kHFNoseSubdetMask)<< kHFNoseSubdetOffset));
}

HFNoseTriggerDetId::HFNoseTriggerDetId(const DetId& gen) {
  if (!gen.null()) {
    if ((gen.det()!=DetId::HGCalTrigger) || gen.subdetId()!=HGCalTriggerSubdetector::HFNoseTrigger) {
      throw cms::Exception("Invalid DetId") << "Cannot initialize HFNoseTriggerDetId from " << std::hex << gen.rawId() << std::dec; 
    }  
  }
  id_ = gen.rawId();
}

HFNoseTriggerDetId& HFNoseTriggerDetId::operator=(const DetId& gen) {
  if (!gen.null()) {
    if ((gen.det()!=DetId::HGCalTrigger) || gen.subdetId()!=HGCalTriggerSubdetector::HFNoseTrigger) {
      throw cms::Exception("Invalid DetId") << "Cannot assign HFNoseTriggerDetId from " << std::hex << gen.rawId() << std::dec; 
    }  
  }
  id_ = gen.rawId();
  return (*this);
}

int HFNoseTriggerDetId::triggerCellX() const {
  // HFNose those need to be redefined
  int N = 12;
  return (3 * (triggerCellV() - N) + 2);

}

int HFNoseTriggerDetId::triggerCellY() const {
  // HFNose those need to be redefined
  int N = 12;
  return (2 * triggerCellU() - (N + triggerCellV()));

}

std::vector<int> HFNoseTriggerDetId::cellU() const {
  // HFNose those need to be redefined
  std::vector<int> uc;
  int nT = 1.;
  if ((triggerCellU() >= HGCalTriggerCell) && 
      (triggerCellV() >= HGCalTriggerCell)) {
    int u0 = nT*triggerCellU(); 
    for (int i=0; i<nT; ++i) {
      for (int j=0; j<nT; ++j) {
	uc.emplace_back(u0+i);
      }
    }
  } else if ((triggerCellU() < HGCalTriggerCell) && 
	     (triggerCellU() <= triggerCellV())) {
    int u0 = nT*triggerCellU(); 
    for (int i=0; i<nT; ++i) {
      for (int j=0; j<nT; ++j) {
	uc.emplace_back(u0+i);
      }
    }
  } else {
    int u0 = nT*(triggerCellU()-1)+1;
    for (int i=0; i<nT; ++i) {
      for (int j=0; j<nT; ++j) {
	uc.emplace_back(u0+j);
      }
      ++u0;
    }
  }
  return uc;
}

std::vector<int> HFNoseTriggerDetId::cellV() const {
  // HFNose those need to be redefined
  std::vector<int> vc;
  int nT = 2.;
  if ((triggerCellU() >= HGCalTriggerCell) && 
      (triggerCellV() >= HGCalTriggerCell)) {
    int v0 = nT*triggerCellV(); 
    for (int i=0; i<nT; ++i) {
      for (int j=0; j<nT; ++j) {
	vc.emplace_back(v0+j);
      }
    }
  } else if ((triggerCellU() < HGCalTriggerCell) && 
	     (triggerCellU() <= triggerCellV())) {
    int v0 = nT*triggerCellV(); 
    for (int i=0; i<nT; ++i) {
      for (int j=0; j<nT; ++j) {
	vc.emplace_back(v0+j);
      }
      ++v0;
    }
  } else {
    int v0 = nT*triggerCellV();
    for (int i=0; i<nT; ++i) {
      for (int j=0; j<nT; ++j) {
	vc.emplace_back(v0+i);
      }
    }
  }
  return vc;
}

std::vector<std::pair<int,int> > HFNoseTriggerDetId::cellUV() const {
  // HFNose those need to be redefined
  std::vector<int> uc = cellU();
  std::vector<int> vc = cellV();
  std::vector<std::pair<int,int> > uv;
  for (unsigned int k=0; k<uc.size(); ++k) {
    uv.emplace_back(std::pair<int,int>(uc[k],vc[k]));
  }
  return uv;
}

std::ostream& operator<<(std::ostream& s,const HFNoseTriggerDetId& id) {
  // HFNose those need to be redefined
  return s << " EE:HSil= " << id.isEE() << ":" << id.isHSilicon()
	   << " type= " << id.type()  << " z= " << id.zside() 
	   << " layer= " << id.layer() 
	   << " wafer(u,v:x,y)= (" << id.waferU() << "," << id.waferV() << ":"
	   << id.waferX() << "," << id.waferY() << ")"
	   << " triggerCell(u,v:x,y)= (" << id.triggerCellU() << "," 
	   << id.triggerCellV() << ":" << id.triggerCellX() << "," 
	   << id.triggerCellY() << ")";
}


