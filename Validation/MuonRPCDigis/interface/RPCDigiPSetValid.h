#ifndef Validation_MuonRPCDigis_RPCDigiPSetValid_h
#define Validation_MuonRPCDigis_RPCDigiPSetValid_h

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPCDigiPSetValid : public DQMEDAnalyzer
{
public:
  RPCDigiPSetValid(const edm::ParameterSet& pset);
  ~RPCDigiPSetValid() override = default;

protected:
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  void bookHistograms(DQMStore::IBooker& booker, edm::Run const& run, edm::EventSetup const& eventSetup) override;

private:
  std::string moduleLabel_;
  bool isFilled_;

  typedef MonitorElement* MEPtr;
  MEPtr hDigiModel1_, hDigiModel2_;

  MEPtr hPropSpeed1_, hPropSpeed2_;
  MEPtr hTimeResol1_, hTimeResol2_;
  MEPtr hTimeJitter1_, hTimeJitter2_;
  MEPtr hEfficiency1_, hEfficiency2_;
  MEPtr hClusterSize1_, hClusterSize2_;
  MEPtr hIRPCTimeResol1_, hIRPCTimeResol2_;
  MEPtr hIRPCTimeJitter1_, hIRPCTimeJitter2_;
};

#endif
