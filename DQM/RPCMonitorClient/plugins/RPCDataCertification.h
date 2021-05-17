#ifndef RPCMonitorClient_RPCDataCertification_H
#define RPCMonitorClient_RPCDataCertification_H

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"

#include <utility>

class RPCDataCertification : public DQMEDHarvester {
public:
  RPCDataCertification(const edm::ParameterSet& pset);
  ~RPCDataCertification() override = default;

protected:
  void beginJob() override;
  void dqmEndLuminosityBlock(DQMStore::IBooker&,
                             DQMStore::IGetter&,
                             edm::LuminosityBlock const&,
                             edm::EventSetup const&) override;      //performed in the endLumi
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override;  //performed in the endJob

private:
  void myBooker(DQMStore::IBooker&);
  void checkFED(edm::EventSetup const&);

  MonitorElement* CertMap_;
  MonitorElement* totalCertFraction;
  constexpr static int kNWheels = 5;
  MonitorElement* certWheelFractions[kNWheels];
  constexpr static int kNDisks = 10;
  MonitorElement* certDiskFractions[kNDisks];
  std::pair<int, int> FEDRange_;
  int numberOfDisks_;
  int NumberOfFeds_;
  bool init_, offlineDQM_;
  double defaultValue_;
};

#endif
