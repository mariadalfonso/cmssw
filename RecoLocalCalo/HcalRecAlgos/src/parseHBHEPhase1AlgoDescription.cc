#include <cassert>

#include "RecoLocalCalo/HcalRecAlgos/interface/parseHBHEPhase1AlgoDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/PulseShapeFitOOTPileupCorrection.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalDeterministicFit.h"

#include "CalibCalorimetry/HcalAlgos/interface/HcalTimeSlew.h"

// Phase 1 HBHE reco algorithm headers
#include "RecoLocalCalo/HcalRecAlgos/interface/SimpleHBHEPhase1Algo.h"

static std::unique_ptr<MahiFit>
parseHBHEMahiDescription(const edm::ParameterSet& conf)
{

  const bool iDynamicPed      = conf.getParameter<bool>   ("dynamicPed");
  const double iTS4Thresh     = conf.getParameter<double> ("ts4Thresh");
  const double chiSqSwitch    = conf.getParameter<double> ("chiSqSwitch");

  const bool iApplyTimeSlew   = conf.getParameter<bool>   ("applyTimeSlew");

  const double iMeanTime      = conf.getParameter<double> ("meanTime");
  const double iTimeSigmaHPD  = conf.getParameter<double> ("timeSigmaHPD");
  const double iTimeSigmaSiPM = conf.getParameter<double> ("timeSigmaSiPM");

  const std::vector<int> iActiveBXs  = conf.getParameter<std::vector<int>> ("activeBXs");
  const int iNMaxItersMin     = conf.getParameter<int>    ("nMaxItersMin");
  const int iNMaxItersNNLS    = conf.getParameter<int>    ("nMaxItersNNLS");
  const double iDeltaChiSqThresh = conf.getParameter<double> ("deltaChiSqThresh");
  const double iNnlsThresh = conf.getParameter<double> ("nnlsThresh");

  std::unique_ptr<MahiFit> corr = std::make_unique<MahiFit>();

  corr->setParameters(iDynamicPed, iTS4Thresh, chiSqSwitch, iApplyTimeSlew, HcalTimeSlew::Medium,
		      iMeanTime, iTimeSigmaHPD, iTimeSigmaSiPM,
		      iActiveBXs, iNMaxItersMin, iNMaxItersNNLS,
		      iDeltaChiSqThresh, iNnlsThresh);

  return corr;
}


static std::unique_ptr<PulseShapeFitOOTPileupCorrection>
parseHBHEMethod2Description(const edm::ParameterSet& conf, const HcalTimeSlew* hcalTimeSlew_delay)
{
    const bool iPedestalConstraint = conf.getParameter<bool>  ("applyPedConstraint");
    const bool iTimeConstraint =     conf.getParameter<bool>  ("applyTimeConstraint");
    const bool iAddPulseJitter =     conf.getParameter<bool>  ("applyPulseJitter");
    const bool iApplyTimeSlew =      conf.getParameter<bool>  ("applyTimeSlew");
    const double iTS4Min =           conf.getParameter<double>("ts4Min");
    const std::vector<double> iTS4Max =           conf.getParameter<std::vector<double>>("ts4Max");
    const double iPulseJitter =      conf.getParameter<double>("pulseJitter");
    const double iTimeMean =         conf.getParameter<double>("meanTime");
    const double iTimeSigHPD =       conf.getParameter<double>("timeSigmaHPD");
    const double iTimeSigSiPM =      conf.getParameter<double>("timeSigmaSiPM");
    const double iPedMean =          conf.getParameter<double>("meanPed");
    const double iTMin =             conf.getParameter<double>("timeMin");
    const double iTMax =             conf.getParameter<double>("timeMax");
    const std::vector<double> its4Chi2 =           conf.getParameter<std::vector<double>>("ts4chi2");
    const int iFitTimes =            conf.getParameter<int>   ("fitTimes");
    if (iTimeConstraint) assert(iTimeSigHPD);
    if (iTimeConstraint) assert(iTimeSigSiPM);

    /////////////////////
    //Sorry for this
    //-C.Madrid
    ////////////////////
    //HcalTimeSlew* hcalTimeSlew_delay = nullptr; 

    std::unique_ptr<PulseShapeFitOOTPileupCorrection> corr =
      std::make_unique<PulseShapeFitOOTPileupCorrection>(hcalTimeSlew_delay);

    corr->setPUParams(iPedestalConstraint, iTimeConstraint, iAddPulseJitter,
                      iApplyTimeSlew, iTS4Min, iTS4Max,
                      iPulseJitter,
		      iTimeMean,iTimeSigHPD, iTimeSigSiPM, iPedMean,
		      iTMin, iTMax, its4Chi2,
                      HcalTimeSlew::Medium, iFitTimes);

    return corr;
}


static std::unique_ptr<HcalDeterministicFit>
parseHBHEMethod3Description(const edm::ParameterSet& conf, const HcalTimeSlew* hcalTimeSlew_delay)
{
    const bool iApplyTimeSlew  =  conf.getParameter<bool>  ("applyTimeSlewM3");
    const float iPedSubThreshold =  conf.getParameter<double>("pedestalUpperLimit");
    const int iTimeSlewParsType  =  conf.getParameter<int>   ("timeSlewParsType");
    const double irespCorrM3 =     conf.getParameter<double>("respCorrM3");
    const std::vector<double>& iTimeSlewPars =
                     conf.getParameter<std::vector<double> >("timeSlewPars");

    PedestalSub pedSubFxn;
    pedSubFxn.init(0, iPedSubThreshold, 0.0);

    //-----------------
    //C. Madrid
    //Need to fix this
    //
    //HcalTimeSlew* hcalTimeSlew_delay = nullptr;

    std::unique_ptr<HcalDeterministicFit> fit = std::make_unique<HcalDeterministicFit>(hcalTimeSlew_delay);
    fit->init( (HcalTimeSlew::ParaSource)iTimeSlewParsType,
	       HcalTimeSlew::Medium, iApplyTimeSlew,
	       pedSubFxn, iTimeSlewPars, irespCorrM3);
    return fit;
}


std::unique_ptr<AbsHBHEPhase1Algo>
parseHBHEPhase1AlgoDescription(const edm::ParameterSet& ps, const HcalTimeSlew* hcalTimeSlew_delay_)
{
    std::unique_ptr<AbsHBHEPhase1Algo> algo;

    const std::string& className = ps.getParameter<std::string>("Class");

    if (className == "SimpleHBHEPhase1Algo")
    {
<<<<<<< HEAD
	std::unique_ptr<MahiFit> mahi;
	std::unique_ptr<PulseShapeFitOOTPileupCorrection> m2;
        std::unique_ptr<HcalDeterministicFit> detFit;

	// only run Mahi OR Method 2 but not both
	if (ps.getParameter<bool>("useMahi") && ps.getParameter<bool>("useM2")) {
          throw cms::Exception("ConfigurationError") <<
            "SimpleHBHEPhase1Algo does not allow both Mahi and Method 2 to be turned on together.";
        }
	if (ps.getParameter<bool>("useMahi"))
	  mahi = parseHBHEMahiDescription(ps);
	if (ps.getParameter<bool>("useM2"))
	  m2 = parseHBHEMethod2Description(ps,hcalTimeSlew_delay_);
	if (ps.getParameter<bool>("useM3"))
	  detFit = parseHBHEMethod3Description(ps,hcalTimeSlew_delay_);

        algo = std::unique_ptr<AbsHBHEPhase1Algo>(
            new SimpleHBHEPhase1Algo(ps.getParameter<int>   ("firstSampleShift"),
                                     ps.getParameter<int>   ("samplesToAdd"),
                                     ps.getParameter<double>("correctionPhaseNS"),
                                     ps.getParameter<double>("tdcTimeShift"),
                                     ps.getParameter<bool>  ("correctForPhaseContainment"),
                                     std::move(m2), std::move(detFit), std::move(mahi))
            );
    }

    return algo;
}
