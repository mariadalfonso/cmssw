#include <algorithm>

#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/SimpleHBHEPhase1Algo.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalCorrectionFunctions.h"

#include "FWCore/Framework/interface/Run.h"

#include "DataFormats/HcalRecHit/interface/HBHERecHitAuxSetter.h"
#include "DataFormats/METReco/interface/HcalPhase1FlagLabels.h"


// Maximum fractional error for calculating Method 0
// pulse containment correction
constexpr float PulseContainmentFractionalError = 0.002f;

SimpleHBHEPhase1Algo::SimpleHBHEPhase1Algo(
    const int firstSampleShift,
    const int samplesToAdd,
    const float phaseNS,
    const float timeShift,
    const bool correctForPhaseContainment,
    const int pulseShapeType, 
    std::unique_ptr<PulseShapeFitOOTPileupCorrection> m2,
    std::unique_ptr<HcalDeterministicFit> detFit,
    std::unique_ptr<DoMahiAlgo> mahi)
    : pulseCorr_(PulseContainmentFractionalError),
      firstSampleShift_(firstSampleShift),
      samplesToAdd_(samplesToAdd),
      phaseNS_(phaseNS),
      timeShift_(timeShift),
      runnum_(0),
      corrFPC_(correctForPhaseContainment),
      pulseShapeType_(pulseShapeType),
      psFitOOTpuCorr_(std::move(m2)),
      hltOOTpuCorr_(std::move(detFit)),
      psFitMAHIOOTpuCorr_(std::move(mahi))
{
  newPulseShapes_ = std::make_unique<NewPulseShapes>();
}

void SimpleHBHEPhase1Algo::beginRun(const edm::Run& r,
                                    const edm::EventSetup& es)
{
    runnum_ = r.run();
    pulseCorr_.beginRun(es);
}

void SimpleHBHEPhase1Algo::endRun()
{
    runnum_ = 0;
    pulseCorr_.endRun();
}

HBHERecHit SimpleHBHEPhase1Algo::reconstruct(const HBHEChannelInfo& info,
                                             const HcalRecoParam* params,
                                             const HcalCalibrations& calibs,
                                             const bool isData)
{

  //temporarily init new pulseshapes
  newPulseShapes_->configurePulseShapes();

    HBHERecHit rh;

    const HcalDetId channelId(info.id());

    // Calculate "Method 0" quantities
    float m0t = 0.f, m0E = 0.f;
    {
        int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int nSamplesToAdd = params ? params->samplesToAdd() : samplesToAdd_;
        const double fc_ampl = info.chargeInWindow(ibeg, ibeg + nSamplesToAdd);
        const bool applyContainment = params ? params->correctForPhaseContainment() : corrFPC_;
        const float phasens = params ? params->correctionPhaseNS() : phaseNS_;
        m0E = m0Energy(info, fc_ampl, applyContainment, phasens, nSamplesToAdd);
        m0E *= hbminusCorrectionFactor(channelId, m0E, isData);
        m0t = m0Time(info, fc_ampl, calibs, nSamplesToAdd);
    }

    // Run "Method 2"
    float m2t = 0.f, m2E = 0.f, chi2 = -1.f;
    bool useTriple = false;
    const PulseShapeFitOOTPileupCorrection* method2 = psFitOOTpuCorr_.get();
    if (method2)
    {

      if(pulseShapeType_==1) {
	// standard pulse shape
	psFitOOTpuCorr_->setPulseShapeTemplate(theHcalPulseShapes_.getShape(info.recoShape()),
					       !info.hasTimeInfo());
      }
      if(pulseShapeType_==2) {
	// effective pulse
	psFitOOTpuCorr_->newSetPulseShapeTemplate(*newPulseShapes_,!info.hasTimeInfo());
	psFitOOTpuCorr_->setDebug(true);
      }

      // "phase1Apply" call below sets m2E, m2t, useTriple, and chi2.
      std::vector<float> tmpResult;
      // These parameters are pased by non-const reference.
      method2->phase1Apply(info, tmpResult, useTriple);
      if (tmpResult.size() != 0) {
	m2t = tmpResult.at(0);
	m2E = tmpResult.at(1);
	chi2 = tmpResult.at(3);
      }
      else {
	m2t=0;
	m2E=0;
	chi2=0;
      }
      m2E *= hbminusCorrectionFactor(channelId, m2E, isData);
    }

    // Run "Method 3"
    float m3t = 0.f, m3E = 0.f;
    const HcalDeterministicFit* method3 = hltOOTpuCorr_.get();
    if (method3)
    {
      if(pulseShapeType_==2) {
	// effective pulse
	hltOOTpuCorr_->configurePulseShapes(*newPulseShapes_);
	hltOOTpuCorr_->setDebug(true);
      }

      // "phase1Apply" sets m3E and m3t (pased by non-const reference)
      std::vector<float> tmpResult;
      method3->phase1Apply(info, tmpResult);
      if (tmpResult.size() != 0) {
	m3E=tmpResult.at(0);
	m3t=tmpResult.at(1);
      }
      else {
	m3E=0;
	m3t=0;
      }
      m3E *= hbminusCorrectionFactor(channelId, m3E, isData);
    }

    // Run "MAHI"
    float m10E = 0.f, chi2_mahi = -1.f;
    DoMahiAlgo* mahi = psFitMAHIOOTpuCorr_.get();
    if (mahi)
    {
      mahi->configurePulseShapes(*newPulseShapes_);

      std::vector<float> tmpResult;
      mahi->phase1Apply(info,tmpResult);
      m10E *= hbminusCorrectionFactor(channelId, m10E, isData);      
      if (tmpResult.size() != 0) {
	m10E=tmpResult.at(0);
	chi2_mahi=tmpResult.at(3);
      }
      else {
	m10E=0;
	chi2_mahi=0;
      }

    }


    // Finally, construct the rechit
    float rhE = m0E;
    float rht = m0t;
    if (method2)
    {
        rhE = m2E;
        rht = m2t;
    }
    else if (method3)
    {
        rhE = m3E;
        rht = m3t;
    }
    else if (mahi)
    {
      rhE = m10E;
      rht = -9999;
      chi2 = chi2_mahi;
    }


    float tdcTime = info.soiRiseTime();
    if (!HcalSpecialTimes::isSpecial(tdcTime))
        tdcTime += timeShift_;
    rh = HBHERecHit(channelId, rhE, rht, tdcTime);
    rh.setRawEnergy(m0E);
    rh.setAuxEnergy(m3E);
    rh.setChiSquared(chi2);

    // Set rechit aux words
    HBHERecHitAuxSetter::setAux(info, &rh);

    // Set some rechit flags (here, for Method 2)
    if (useTriple)
       rh.setFlagField(1, HcalPhase1FlagLabels::HBHEPulseFitBit);

    return rh;
}

float SimpleHBHEPhase1Algo::hbminusCorrectionFactor(const HcalDetId& cell,
                                                    const float energy,
                                                    const bool isRealData) const
{
    float corr = 1.f;
    if (isRealData && runnum_ > 0)
        if (cell.subdet() == HcalBarrel)
        {
            const int ieta = cell.ieta();
            const int iphi = cell.iphi();
            corr = hbminus_special_ecorr(ieta, iphi, energy, runnum_);
        }
    return corr;
}

float SimpleHBHEPhase1Algo::m0Energy(const HBHEChannelInfo& info,
                                     const double fc_ampl,
                                     const bool applyContainmentCorrection,
                                     const double phaseNs,
                                     const int nSamplesToAdd)
{
    int ibeg = static_cast<int>(info.soi()) + firstSampleShift_;
    if (ibeg < 0)
        ibeg = 0;
    double e = info.energyInWindow(ibeg, ibeg + nSamplesToAdd);

    // Pulse containment correction
    {    
        double corrFactor = 1.0;
        if (applyContainmentCorrection)
            corrFactor = pulseCorr_.get(info.id(), nSamplesToAdd, phaseNs)->getCorrection(fc_ampl);
        e *= corrFactor;
    }

    return e;
}

float SimpleHBHEPhase1Algo::m0Time(const HBHEChannelInfo& info,
                                   const double fc_ampl,
                                   const HcalCalibrations& calibs,
                                   const int nSamplesToExamine) const
{
    float time = -9999.f; // historic value

    const unsigned nSamples = info.nSamples();
    if (nSamples > 2U)
    {
        const int soi = info.soi();
        int ibeg = soi + firstSampleShift_;
        if (ibeg < 0)
            ibeg = 0;
        const int iend = ibeg + nSamplesToExamine;
        unsigned maxI = info.peakEnergyTS(ibeg, iend);
        if (maxI < HBHEChannelInfo::MAXSAMPLES)
        {
            if (!maxI)
                maxI = 1U;
            else if (maxI >= nSamples - 1U)
                maxI = nSamples - 2U;

            // The remaining code in this scope emulates
            // the historic algorithm
            float t0 = info.tsEnergy(maxI - 1U);
            float maxA = info.tsEnergy(maxI);
            float t2 = info.tsEnergy(maxI + 1U);

            // Handle negative excursions by moving "zero"
            float minA = t0;
            if (maxA < minA) minA = maxA;
            if (t2 < minA)   minA=t2;
            if (minA < 0.f) { maxA-=minA; t0-=minA; t2-=minA; }
            float wpksamp = (t0 + maxA + t2);
            if (wpksamp) wpksamp = (maxA + 2.f*t2) / wpksamp;
            time = (maxI - soi)*25.f + timeshift_ns_hbheho(wpksamp);

            // Legacy QIE8 timing correction
            time -= HcalTimeSlew::delay(std::max(1.0, fc_ampl),
                                        HcalTimeSlew::Medium);
            // Time calibration
            time -= calibs.timecorr();
        }
    }
    return time;
}
