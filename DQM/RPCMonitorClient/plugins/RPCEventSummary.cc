/*  \author Anna Cimmino*/
#include <DQM/RPCMonitorClient/plugins/RPCEventSummary.h>

//CondFormats
#include "CondFormats/RunInfo/interface/RunInfo.h"
#include "CondFormats/DataRecord/interface/RunSummaryRcd.h"
// Framework
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include <fmt/format.h>

RPCEventSummary::RPCEventSummary(const edm::ParameterSet& ps) {
  edm::LogVerbatim("rpceventsummary") << "[RPCEventSummary]: Constructor";

  enableReportSummary_ = ps.getUntrackedParameter<bool>("EnableSummaryReport", true);
  prescaleFactor_ = ps.getUntrackedParameter<int>("PrescaleFactor", 1);
  eventInfoPath_ = ps.getUntrackedParameter<std::string>("EventInfoPath", "RPC/EventInfo");

  std::string subsystemFolder = ps.getUntrackedParameter<std::string>("RPCFolder", "RPC");
  std::string recHitTypeFolder = ps.getUntrackedParameter<std::string>("RecHitTypeFolder", "AllHits");
  std::string summaryFolder = ps.getUntrackedParameter<std::string>("SummaryFolder", "SummaryHistograms");

  globalFolder_ = subsystemFolder + "/" + recHitTypeFolder + "/" + summaryFolder;
  prefixFolder_ = subsystemFolder + "/" + recHitTypeFolder;

  minimumEvents_ = ps.getUntrackedParameter<int>("MinimumRPCEvents", 10000);
  numberDisk_ = ps.getUntrackedParameter<int>("NumberOfEndcapDisks", 4);
  doEndcapCertification_ = ps.getUntrackedParameter<bool>("EnableEndcapSummary", false);

  FEDRange_.first = ps.getUntrackedParameter<unsigned int>("MinimumRPCFEDId", 790);
  FEDRange_.second = ps.getUntrackedParameter<unsigned int>("MaximumRPCFEDId", 792);

  NumberOfFeds_ = FEDRange_.second - FEDRange_.first + 1;

  offlineDQM_ = ps.getUntrackedParameter<bool>("OfflineDQM", true);
}

RPCEventSummary::~RPCEventSummary() { edm::LogVerbatim("rpceventsummary") << "[RPCEventSummary]: Destructor "; }

void RPCEventSummary::beginJob() {
  edm::LogVerbatim("rpceventsummary") << "[RPCEventSummary]: Begin job ";
  init_ = false;
}

void RPCEventSummary::dqmEndLuminosityBlock(DQMStore::IBooker& ibooker,
                                            DQMStore::IGetter& igetter,
                                            edm::LuminosityBlock const& lb,
                                            edm::EventSetup const& setup) {
  edm::LogVerbatim("rpceventsummary") << "[RPCEventSummary]: Begin run";

  if (!init_) {
    lumiCounter_ = prescaleFactor_;

    edm::eventsetup::EventSetupRecordKey recordKey(
        edm::eventsetup::EventSetupRecordKey::TypeTag::findType("RunInfoRcd"));

    int defaultValue = 1;
    isIn_ = true;

    if (auto runInfoRec = setup.tryToGet<RunInfoRcd>()) {
      defaultValue = -1;
      //get fed summary information
      edm::ESHandle<RunInfo> sumFED;
      runInfoRec->get(sumFED);
      std::vector<int> FedsInIds = sumFED->m_fed_in;
      unsigned int f = 0;
      isIn_ = false;
      while (!isIn_ && f < FedsInIds.size()) {
        int fedID = FedsInIds[f];
        //make sure fed id is in allowed range
        if (fedID >= FEDRange_.first && fedID <= FEDRange_.second) {
          defaultValue = 1;
          isIn_ = true;
        }
        f++;
      }
    }

    MonitorElement* me;
    ibooker.setCurrentFolder(eventInfoPath_);

    //a global summary float [0,1] providing a global summary of the status
    //and showing the goodness of the data taken by the the sub-system
    std::string histoName = "reportSummary";
    me = nullptr;
    me = ibooker.bookFloat(histoName);
    me->Fill(defaultValue);

    //TH2F ME providing a mapof values[0-1] to show if problems are localized or distributed
    me = nullptr;
    me = ibooker.book2D("reportSummaryMap", "RPC Report Summary Map", 15, -7.5, 7.5, 12, 0.5, 12.5);

    //customize the 2d histo
    for (int i = 1; i <= 15; i++) {
      if (i < 13) {
        me->setBinLabel(i, fmt::format("Sec{}", i), 2);
      }

      std::string binLabel;
      if (i < 5)
        binLabel = fmt::format("Disk{}", i-5);
      else if (i > 11)
        binLabel = fmt::format("Disk{}", i-11);
      else if (i != 11 and i != 5)
        binLabel = fmt::format("Wheel{}", i-8);

      me->setBinLabel(i, binLabel, 1);
    }

    //fill the histo with "1" --- just for the moment
    for (int i = 1; i <= 15; i++) {
      for (int j = 1; j <= 12; j++) {
        if (i == 5 || i == 11 || (j > 6 && (i < 6 || i > 10)))
          me->setBinContent(i, j, -1);  //bins that not correspond to subdetector parts
        else
          me->setBinContent(i, j, defaultValue);
      }
    }

    if (numberDisk_ < 4)
      for (int j = 1; j <= 12; j++) {
        me->setBinContent(1, j, -1);  //bins that not correspond to subdetector parts
        me->setBinContent(15, j, -1);
      }

    //the reportSummaryContents folder containins a collection of ME floats [0-1] (order of 5-10)
    // which describe the behavior of the respective subsystem sub-components.
    ibooker.setCurrentFolder(eventInfoPath_ + "/reportSummaryContents");

    std::vector<std::string> segmentNames;
    for (int i = -2; i <= 2; i++) {
      segmentNames.push_back(fmt::format("RPC_Wheel{}", i));
    }

    for (int i = -numberDisk_; i <= numberDisk_; i++) {
      if (i == 0)
        continue;
      segmentNames.push_back(fmt::format("RPC_Disk{}", i));
    }

    for (unsigned int i = 0; i < segmentNames.size(); i++) {
      me = ibooker.bookFloat(segmentNames[i]);
      me->Fill(defaultValue);
    }

    lumiCounter_ = prescaleFactor_;
    init_ = true;
  }

  if (isIn_ && !offlineDQM_ && lumiCounter_ % prescaleFactor_ == 0) {
    this->clientOperation(igetter);
  }

  lumiCounter_++;
}

void RPCEventSummary::dqmEndJob(DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter) {
  if (isIn_) {
    this->clientOperation(igetter);
  }
}

void RPCEventSummary::clientOperation(DQMStore::IGetter& igetter) {
  float rpcevents = minimumEvents_;
  MonitorElement* meRPCEvents = igetter.get(prefixFolder_ + "/RPCEvents");

  if (meRPCEvents) {
    rpcevents = meRPCEvents->getBinContent(1);
  }

  if (rpcevents < minimumEvents_)
    return;

  MonitorElement* reportMe = igetter.get(eventInfoPath_ + "/reportSummaryMap");

  //BARREL
  float barrelFactor = 0;
  float endcapFactor = 0;

  //Fill repor summary
  float rpcFactor = barrelFactor;
  if (doEndcapCertification_) {
    rpcFactor = (barrelFactor + endcapFactor) / 2;
  }

  MonitorElement* globalMe = igetter.get(eventInfoPath_ + "/reportSummary");
  if (globalMe)
    globalMe->Fill(rpcFactor);
}
