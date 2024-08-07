#include "DQM/EcalMonitorClient/interface/LedClient.h"

#include "DataFormats/EcalDetId/interface/EcalPnDiodeDetId.h"

#include "CondFormats/EcalObjects/interface/EcalDQMStatusHelper.h"

#include "DQM/EcalCommon/interface/EcalDQMCommonUtils.h"
#include "DQM/EcalCommon/interface/MESetMulti.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include <cmath>
#include <fstream>
namespace ecaldqm {
  LedClient::LedClient()
      : DQWorkerClient(),
        wlToME_(),
        minChannelEntries_(0),
        expectedAmplitude_(0),
        toleranceAmplitude_(0.),
        toleranceAmpRMSRatio_(0.),
        expectedTiming_(0),
        toleranceTiming_(0.),
        toleranceTimRMS_(0.),
        expectedPNAmplitude_(0),
        tolerancePNAmp_(0.),
        tolerancePNRMSRatio_(0.),
        forwardFactor_(0.) {}

  void LedClient::setParams(edm::ParameterSet const& _params) {
    minChannelEntries_ = _params.getUntrackedParameter<int>("minChannelEntries");
    toleranceAmplitude_ = _params.getUntrackedParameter<double>("toleranceAmplitude");
    toleranceAmpRMSRatio_ = _params.getUntrackedParameter<double>("toleranceAmpRMSRatio");
    toleranceTiming_ = _params.getUntrackedParameter<double>("toleranceTiming");
    toleranceTimRMS_ = _params.getUntrackedParameter<double>("toleranceTimRMS");
    tolerancePNAmp_ = _params.getUntrackedParameter<double>("tolerancePNAmp");
    tolerancePNRMSRatio_ = _params.getUntrackedParameter<double>("tolerancePNRMSRatio");
    forwardFactor_ = _params.getUntrackedParameter<double>("forwardFactor");

    std::vector<int> ledWavelengths(_params.getUntrackedParameter<std::vector<int> >("ledWavelengths"));

    // wavelengths are not necessarily ordered
    // create a map wl -> MESet index
    // using Amplitude here but any multi-wavelength plot is fine

    MESet::PathReplacements repl;

    MESetMulti const& amplitude(static_cast<MESetMulti const&>(sources_.at("Amplitude")));
    unsigned nWL(ledWavelengths.size());
    for (unsigned iWL(0); iWL != nWL; ++iWL) {
      int wl(ledWavelengths[iWL]);
      if (wl != 1 && wl != 2)
        throw cms::Exception("InvalidConfiguration") << "Led Wavelength";
      repl["wl"] = std::to_string(wl);
      wlToME_[wl] = amplitude.getIndex(repl);
    }

    expectedAmplitude_.resize(nWL);
    expectedTiming_.resize(nWL);
    expectedPNAmplitude_.resize(nWL);

    std::vector<double> inExpectedAmplitude(_params.getUntrackedParameter<std::vector<double> >("expectedAmplitude"));
    std::vector<double> inExpectedTiming(_params.getUntrackedParameter<std::vector<double> >("expectedTiming"));
    std::vector<double> inExpectedPNAmplitude(
        _params.getUntrackedParameter<std::vector<double> >("expectedPNAmplitude"));

    for (std::map<int, unsigned>::iterator wlItr(wlToME_.begin()); wlItr != wlToME_.end(); ++wlItr) {
      unsigned iME(wlItr->second);
      int iWL(wlItr->first - 1);
      expectedAmplitude_[iME] = inExpectedAmplitude[iWL];
      expectedTiming_[iME] = inExpectedTiming[iWL];
      expectedPNAmplitude_[iME] = inExpectedPNAmplitude[iWL];
    }

    //Get the list of known problematic Supercrystal ids and store them in the vector SClist_
    std::string SClistpath = edm::FileInPath("DQM/EcalMonitorClient/data/LedTowers/SClist.dat").fullPath();
    std::ifstream infile;
    infile.open((SClistpath).c_str());
    uint32_t detid;
    int ix, iy, iz;
    while (!infile.eof()) {
      infile >> ix >> iy >> iz >> detid;
      SClist_.push_back(detid);
    }

    qualitySummaries_.insert("Quality");
    qualitySummaries_.insert("QualitySummary");
    qualitySummaries_.insert("PNQualitySummary");
  }

  void LedClient::producePlots(ProcessType) {
    uint32_t mask(1 << EcalDQMStatusHelper::LED_MEAN_ERROR | 1 << EcalDQMStatusHelper::LED_RMS_ERROR |
                  1 << EcalDQMStatusHelper::LED_TIMING_MEAN_ERROR | 1 << EcalDQMStatusHelper::LED_TIMING_RMS_ERROR);

    MESetMulti& meQuality(static_cast<MESetMulti&>(MEs_.at("Quality")));
    MESetMulti& meQualitySummary(static_cast<MESetMulti&>(MEs_.at("QualitySummary")));
    MESetMulti& meAmplitudeMean(static_cast<MESetMulti&>(MEs_.at("AmplitudeMean")));
    MESetMulti& meAmplitudeRMS(static_cast<MESetMulti&>(MEs_.at("AmplitudeRMS")));
    MESetMulti& meTimingMean(static_cast<MESetMulti&>(MEs_.at("TimingMean")));
    MESetMulti& meTimingRMSMap(static_cast<MESetMulti&>(MEs_.at("TimingRMSMap")));
    MESetMulti& mePNQualitySummary(static_cast<MESetMulti&>(MEs_.at("PNQualitySummary")));

    MESetMulti const& sAmplitude(static_cast<MESetMulti const&>(sources_.at("Amplitude")));
    MESetMulti const& sTiming(static_cast<MESetMulti const&>(sources_.at("Timing")));
    MESetMulti const& sPNAmplitude(static_cast<MESetMulti const&>(sources_.at("PNAmplitude")));
    MESet const& sCalibStatus(static_cast<MESet const&>(sources_.at("CalibStatus")));

    for (std::map<int, unsigned>::iterator wlItr(wlToME_.begin()); wlItr != wlToME_.end(); ++wlItr) {
      meQuality.use(wlItr->second);
      meQualitySummary.use(wlItr->second);
      meAmplitudeMean.use(wlItr->second);
      meAmplitudeRMS.use(wlItr->second);
      meTimingMean.use(wlItr->second);
      meTimingRMSMap.use(wlItr->second);
      mePNQualitySummary.use(wlItr->second);

      sAmplitude.use(wlItr->second);
      sTiming.use(wlItr->second);
      sPNAmplitude.use(wlItr->second);

      MESet::iterator qEnd(meQuality.end(GetElectronicsMap()));

      MESet::const_iterator tItr(GetElectronicsMap(), sTiming);
      MESet::const_iterator aItr(GetElectronicsMap(), sAmplitude);

      int wl(wlItr->first + 3);
      bool enabled(wl < 0 ? false : sCalibStatus.getBinContent(getEcalDQMSetupObjects(), wl) > 0 ? true : false);
      for (MESet::iterator qItr(meQuality.beginChannel(GetElectronicsMap())); qItr != qEnd;
           qItr.toNextChannel(GetElectronicsMap())) {
        DetId id(qItr->getId());

        bool doMask(meQuality.maskMatches(id, mask, statusManager_, GetTrigTowerMap()));

        aItr = qItr;

        float aEntries(aItr->getBinEntries());

        if (aEntries < minChannelEntries_) {
          qItr->setBinContent(enabled ? (doMask ? kMUnknown : kUnknown) : kMUnknown);
          continue;
        }

        float aMean(aItr->getBinContent());
        float aRms(aItr->getBinError() * sqrt(aEntries));

        meAmplitudeMean.fill(getEcalDQMSetupObjects(), id, aMean);
        meAmplitudeRMS.setBinContent(getEcalDQMSetupObjects(), id, aRms);

        tItr = qItr;

        float tEntries(tItr->getBinEntries());

        if (tEntries < minChannelEntries_)
          continue;

        float tMean(tItr->getBinContent());
        float tRms(tItr->getBinError() * sqrt(tEntries));

        meTimingMean.fill(getEcalDQMSetupObjects(), id, tMean);
        meTimingRMSMap.setBinContent(getEcalDQMSetupObjects(), id, tRms);

        //Temporarily disabling all cuts on LED Quality plot.
        qItr->setBinContent(doMask ? kMGood : kGood);

        /*
        float intensity(aMean / expectedAmplitude_[wlItr->second]);
        if (isForward(id))
          intensity /= forwardFactor_;

        float aRmsThr(sqrt(pow(aMean * toleranceAmpRMSRatio_, 2) + pow(3., 2)));

        EcalScDetId scid = EEDetId(id).sc();  //Get the Endcap SC id for the given crystal id.

        //For the known bad Supercrystals in the SClist, bad quality flag is only set based on the amplitude RMS
        //and everything else is ignored.
        if (std::find(SClist_.begin(), SClist_.end(), int(scid)) != SClist_.end()) {
          if (aRms > aRmsThr)
            qItr->setBinContent(doMask ? kMBad : kBad);
          else
            qItr->setBinContent(doMask ? kMGood : kGood);
        } else {
          if (intensity < toleranceAmplitude_ || aRms > aRmsThr ||
              std::abs(tMean - expectedTiming_[wlItr->second]) > toleranceTiming_ || tRms > toleranceTimRMS_)
            qItr->setBinContent(doMask ? kMBad : kBad);
          else
            qItr->setBinContent(doMask ? kMGood : kGood);
        }*/
      }

      towerAverage_(meQualitySummary, meQuality, 0.2);

      for (unsigned iDCC(0); iDCC < nDCC; ++iDCC) {
        if (memDCCIndex(iDCC + 1) == unsigned(-1))
          continue;
        if (iDCC >= kEBmLow && iDCC <= kEBpHigh)
          continue;

        for (unsigned iPN(0); iPN < 10; ++iPN) {
          EcalPnDiodeDetId id(EcalEndcap, iDCC + 1, iPN + 1);

          bool doMask(mePNQualitySummary.maskMatches(id, mask, statusManager_, GetTrigTowerMap()));

          float pEntries(sPNAmplitude.getBinEntries(getEcalDQMSetupObjects(), id));

          if (pEntries < minChannelEntries_) {
            mePNQualitySummary.setBinContent(getEcalDQMSetupObjects(), id, doMask ? kMUnknown : kUnknown);
            continue;
          }

          float pMean(sPNAmplitude.getBinContent(getEcalDQMSetupObjects(), id));
          float pRms(sPNAmplitude.getBinError(getEcalDQMSetupObjects(), id) * sqrt(pEntries));
          float intensity(pMean / expectedPNAmplitude_[wlItr->second]);

          if (intensity < tolerancePNAmp_ || pRms > pMean * tolerancePNRMSRatio_)
            mePNQualitySummary.setBinContent(getEcalDQMSetupObjects(), id, doMask ? kMBad : kBad);
          else
            mePNQualitySummary.setBinContent(getEcalDQMSetupObjects(), id, doMask ? kMGood : kGood);
        }
      }
    }
  }

  DEFINE_ECALDQM_WORKER(LedClient);
}  // namespace ecaldqm
