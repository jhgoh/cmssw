#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "FWCore/Framework/interface/Event.h"
//#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
//#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

class RPCEfficiencyTagProbe : public DQMEDAnalyzer
{
public:
  RPCEfficiencyTagProbe(const edm::ParameterSet&);
  ~RPCEfficiencyTagProbe() override = default;

protected:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void dqmBeginRun(const edm::Run&, const edm::EventSetup&) override;
  void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;

  constexpr static double muonMass_ = 0.1056583715;

private:
  // Objects to build tag and probe pair
  const edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  const edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  const edm::EDGetTokenT<reco::MuonCollection> muonToken_;
  const edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;
  const edm::EDGetTokenT<RPCRecHitCollection> rpcHitToken_;

  const double minMuonPt_, maxMuonAbsEta_, maxMuonRelIso_;
  const double minTrackPt_, maxTrackAbsEta_;
  const bool doCheckSign_;
  const double minMass_, maxMass_;
  const double minDR_;
  const std::vector<std::string> triggerPaths_, triggerModules_;

  HLTConfigProvider hltConfig_;

private:
  // Monitor elements
  typedef MonitorElement* MEP;
  MEP h_tpMass_;
  MEP h_tagPt_, h_probePt_;

  std::map<long, int> rawIdToBin_;
  MEP h_detIdAll_, h_detIdPass_;
};

