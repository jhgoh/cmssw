#ifndef RPCOccupancyTest_H
#define RPCOccupancyTest_H

#include "DQM/RPCMonitorClient/interface/RPCClient.h"

class RPCOccupancyTest : public RPCClient {
public:
  RPCOccupancyTest(const edm::ParameterSet &ps);
  ~RPCOccupancyTest() override = default;

  void clientOperation() override;
  void getMonitorElements(std::vector<MonitorElement *> &, std::vector<RPCDetId> &, std::string &) override;
  void beginJob(std::string &) override;
  void myBooker(DQMStore::IBooker &) override;

protected:
  // void OccupancyDist();
  void fillGlobalME(RPCDetId &, MonitorElement *);

private:
  std::string globalFolder_, prefixDir_;
  std::vector<MonitorElement *> myOccupancyMe_;
  std::vector<RPCDetId> myDetIds_;

  float totalActive_, totalStrips_;

  MonitorElement *Active_Fraction;  // Fraction of channels with data
  MonitorElement *Active_Dead;
};

#endif
