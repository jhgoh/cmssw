#ifndef Validation_MuonRPCDigis_RPCSimDigiValid_h
#define Validation_MuonRPCDigis_RPCSimDigiValid_h

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

class RPCSimDigiValid : public DQMEDAnalyzer
{
public:
  RPCSimDigiValid(const edm::ParameterSet& pset);
  ~RPCSimDigiValid() override = default;

protected:
  void analyze(const edm::Event& event, const edm::EventSetup& eventSetup) override;
  void bookHistograms(DQMStore::IBooker& booker, edm::Run const& run, edm::EventSetup const& eventSetup) override;

private:
  constexpr static double speedOfLight_ = 29.9792458; // spped of light in [cm/ns]

  edm::EDGetTokenT<edm::PSimHitContainer> simHitsToken_;
  edm::EDGetTokenT<RPCDigiCollection> rpcDigisToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;


  typedef MonitorElement* MEPtr;

  MEPtr hNEvents_;
  MEPtr hNPileup_;

  MEPtr hSimHitPid_;
  MEPtr hSimHitTime_;
  MEPtr hSimHitTime_Mu_, hSimHitTime_El_, hSimHitTime_Others_;
  MEPtr hSimHitELoss_;
  MEPtr hSimHitELoss_Mu_, hSimHitELoss_El_, hSimHitELoss_Others_;

  MEPtr hSimHitELoss_Barrel_;
  MEPtr hSimHitR_Barrel_;
  MEPtr hSimHitZ_Barrel_;
  MEPtr hSimHitZ_RB1in_, hSimHitZ_RB1out_, hSimHitZ_RB2in_, hSimHitZ_RB2out_, hSimHitZ_RB3_, hSimHitZ_RB4_;
  MEPtr hSimHitZ_RB1in_Mu_, hSimHitZ_RB1out_Mu_, hSimHitZ_RB2in_Mu_, hSimHitZ_RB2out_Mu_, hSimHitZ_RB3_Mu_, hSimHitZ_RB4_Mu_;

  MEPtr hSimHitELoss_EndcapP_;
  MEPtr hSimHitR_EndcapP_;
  MEPtr hSimHitR_REP1_, hSimHitR_REP2_, hSimHitR_REP3_, hSimHitR_REP4_;
  MEPtr hSimHitR_REP1_Mu_, hSimHitR_REP2_Mu_, hSimHitR_REP3_Mu_, hSimHitR_REP4_Mu_;
  MEPtr hSimHitZ_EndcapP_;

  MEPtr hSimHitELoss_EndcapM_;
  MEPtr hSimHitR_EndcapM_;
  MEPtr hSimHitR_REM1_, hSimHitR_REM2_, hSimHitR_REM3_, hSimHitR_REM4_;
  MEPtr hSimHitR_REM1_Mu_, hSimHitR_REM2_Mu_, hSimHitR_REM3_Mu_, hSimHitR_REM4_Mu_;
  MEPtr hSimHitZ_EndcapM_;
};

#endif
