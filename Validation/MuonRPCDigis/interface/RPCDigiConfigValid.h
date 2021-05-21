#ifndef Validation_MuonRPCDigis_RPCDigiConfigValid_h
#define Validation_MuonRPCDigis_RPCDigiConfigValid_h

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

class RPCDigiConfigValid : public DQMEDAnalyzer
{
public:
  RPCDigiConfigValid(const edm::ParameterSet& pset);
  ~RPCDigiConfigValid() override = default;

protected:
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  void bookHistograms(DQMStore::IBooker& booker, edm::Run const& run, edm::EventSetup const& eventSetup) override;

private:
  std::string moduleLabel_;
  bool isFilled_;

  typedef MonitorElement* MEPtr;
  MEPtr hDigiModel1_, hDigiModel2_;
  MEPtr hNoise_, hSignal_, hDoBkgNoise_;

  MEPtr hPropSpeed1_, hPropSpeed2_;
  MEPtr hTimeResol1_, hTimeResol2_;
  MEPtr hTimeJitter1_, hTimeJitter2_;
  MEPtr hEfficiency1_, hEfficiency2_;
  MEPtr hClusterSize1_, hClusterSize2_;
  MEPtr hDigiElectron1_, hDigiElectron2_;
  MEPtr hLinkGateWidth1_, hLinkGateWidth2_;
  MEPtr hIRPCTimeResol1_, hIRPCTimeResol2_;
  MEPtr hIRPCTimeJitter1_, hIRPCTimeJitter2_;
};

#endif
