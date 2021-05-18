/*  \author Anna Cimmino*/
#include "DQM/RPCMonitorClient/interface/RPCOccupancyTest.h"
#include "DQM/RPCMonitorDigi/interface/utils.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"

#include <fmt/format.h>

RPCOccupancyTest::RPCOccupancyTest(const edm::ParameterSet& ps) {
  edm::LogVerbatim("rpceventsummary") << "[RPCOccupancyTest]: Constructor";
}

void RPCOccupancyTest::beginJob(std::string& workingFolder) {
  edm::LogVerbatim("rpceventsummary") << "[RPCOccupancyTest]: Begin job ";
  globalFolder_ = workingFolder;

  totalStrips_ = 0.;
  totalActive_ = 0.;
}

void RPCOccupancyTest::getMonitorElements(std::vector<MonitorElement*>& meVector,
                                          std::vector<RPCDetId>& detIdVector,
                                          std::string& clientHistoName) {
  //Get NumberOfDigi ME for each roll
  for (unsigned int i = 0; i < meVector.size(); i++) {
    const std::string meName = meVector[i]->getName();

    if (meName.find(clientHistoName) != std::string::npos) {
      myOccupancyMe_.push_back(meVector[i]);
      myDetIds_.push_back(detIdVector[i]);
    }
  }
}

void RPCOccupancyTest::clientOperation() {
  edm::LogVerbatim("rpceventsummary") << "[RPCOccupancyTest]: Client Operation";

  //Loop on MEs
  for (unsigned int i = 0; i < myOccupancyMe_.size(); i++) {
    this->fillGlobalME(myDetIds_[i], myOccupancyMe_[i]);
  }  //End loop on MEs

  //Active Channels
  if (Active_Fraction && totalStrips_ != 0.) {
    Active_Fraction->setBinContent(1, (totalActive_ / totalStrips_));
  }
  if (Active_Dead) {
    Active_Dead->setBinContent(1, totalActive_);
    Active_Dead->setBinContent(2, (totalStrips_ - totalActive_));
  }
}

void RPCOccupancyTest::myBooker(DQMStore::IBooker& ibooker) {
  ibooker.setCurrentFolder(globalFolder_);

  Active_Fraction = ibooker.book1D("RPC_Active_Channel_Fractions", "RPC_Active_Channel_Fractions", 1, 0.5, 1.5);
  Active_Fraction->setBinLabel(1, "Active Fraction", 1);

  Active_Dead = ibooker.book1D("RPC_Active_Inactive_Strips", "RPC_Active_Inactive_Strips", 2, 0.5, 2.5);
  Active_Dead->setBinLabel(1, "Active Strips", 1);
  Active_Dead->setBinLabel(2, "Inactive Strips", 1);
}

void RPCOccupancyTest::fillGlobalME(RPCDetId& detId, MonitorElement* myMe) {
  if (!myMe)
    return;

  const int stripInRoll = myMe->getNbinsX();
  totalStrips_ += stripInRoll;

  for (int strip = 1; strip <= stripInRoll; strip++) {
    const float stripEntries = myMe->getBinContent(strip);
    if (stripEntries > 0) {
      ++totalActive_;
    }
  }
}
