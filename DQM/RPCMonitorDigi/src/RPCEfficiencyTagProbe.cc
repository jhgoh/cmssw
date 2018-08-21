#include "DQM/RPCMonitorDigi/interface/RPCEfficiencyTagProbe.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"

#include "DataFormats/Math/interface/deltaR.h"

#include <vector>
#include <string>
#include <memory>
#include <cmath>

RPCEfficiencyTagProbe::RPCEfficiencyTagProbe(const edm::ParameterSet& pset):
  triggerResultsToken_(consumes<edm::TriggerResults>(pset.getParameter<edm::InputTag>("triggerResults"))),
  pvToken_(consumes<reco::VertexCollection>(pset.getParameter<edm::InputTag>("vertex"))),
  muonToken_(consumes<reco::MuonCollection>(pset.getParameter<edm::InputTag>("muons"))),
  triggerEventToken_(consumes<trigger::TriggerEvent>(pset.getParameter<edm::InputTag>("triggerObjects"))),
  rpcHitToken_(consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcRecHits"))),
  minMuonPt_(pset.getParameter<double>("minMuonPt")),
  maxMuonAbsEta_(pset.getParameter<double>("maxMuonAbsEta")),
  maxMuonRelIso_(pset.getParameter<double>("maxMuonRelIso")),
  minTrackPt_(pset.getParameter<double>("minTrackPt")),
  maxTrackAbsEta_(pset.getParameter<double>("maxTrackAbsEta")),
  doCheckSign_(pset.getParameter<bool>("doCheckSign")),
  minMass_(pset.getParameter<double>("minMass")),
  maxMass_(pset.getParameter<double>("maxMass")),
  minDR_(pset.getParameter<double>("minDR")),
  triggerPaths_(pset.getParameter<std::vector<std::string>>("triggerPaths")),
  triggerModules_(pset.getParameter<std::vector<std::string>>("triggerModules"))
{
}

void RPCEfficiencyTagProbe::dqmBeginRun(const edm::Run& run, const edm::EventSetup& eventSetup)
{
  bool changed = true;
  hltConfig_.init(run, eventSetup, "HLT", changed);
}

void RPCEfficiencyTagProbe::bookHistograms(DQMStore::IBooker& booker, edm::Run const&, edm::EventSetup const& eventSetup)
{
  booker.setCurrentFolder("RPC/RPCEfficiency/TagAndProbe");

  h_tpMass_ = booker.book1D("tpMass", "Mass of tag-probe pair;Invariant mass (GeV);Events", 100, minMass_, maxMass_);
  h_tagPt_ = booker.book1D("tagPt", "Transverse momentum of tag muon pt;Transverse momentum (GeV);Events", 100, 0, 100);
  h_probePt_ = booker.book1D("probePt", "Transverse momentum of probe muon pt;Transverse momentum (GeV);Events", 100, 0, 100);

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);
  std::vector<std::string> rollNames;
  rollNames.reserve(rpcGeom->rolls().size());
  for ( auto roll : rpcGeom->rolls() ) {
    const int rawId = roll->geographicalId().rawId();
    RPCDetId rpcId(rawId);

    rawIdToBin_[rawId] = rollNames.size()+1; // bin number starts from 1
    rollNames.push_back(RPCGeomServ(rpcId).name());
  }
  const size_t nRolls = rollNames.size();

  // Book hit counts vs rawId histograms. We do not add bin labels to save memory.
  // We can recover the name of rolls and fill appropriate plots and bins during the harvesting step.
  h_detIdAll_ = booker.book1D("detIdAll", "RPC DetId with valid extrapolation;DetId;Events", nRolls, 1, nRolls+1);
  h_detIdPass_ = booker.book1D("detIdPass", "RPC DetId with matched RecHit;DetId;Events", nRolls, 1, nRolls+1);

  // Set bin labels
  auto axisDetIdAll = h_detIdAll_->getTH1()->GetXaxis();
  auto axisDetIdPass = h_detIdPass_->getTH1()->GetXaxis();
  for ( size_t i=0; i<nRolls; ++i ) {
    axisDetIdAll->SetBinLabel(i+1, rollNames[i].c_str());
    axisDetIdPass->SetBinLabel(i+1, rollNames[i].c_str());
  }
}

void RPCEfficiencyTagProbe::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  double mass = -1;

  edm::Handle<edm::TriggerResults> triggerResultsHandle;
  event.getByToken(triggerResultsToken_, triggerResultsHandle);

  edm::Handle<reco::VertexCollection> pvHandle;
  event.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::MuonCollection> muonHandle;
  event.getByToken(muonToken_, muonHandle);

  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  event.getByToken(triggerEventToken_, triggerEventHandle);

  reco::MuonRef tagRef;
  reco::MuonRef probeRef;
  const reco::Vertex* goodPV = nullptr;
  for ( auto& pv : *pvHandle ) {
    if ( !pv.isFake() and pv.ndof() >= 4.0 and
         std::abs(pv.z()) <= 24.0 and std::abs(pv.position().Rho()) <= 2.0 ) {
      goodPV = &pv;
      break;
    }
  }
  if ( goodPV == nullptr ) return;

  // Collect interested trigger objects
  auto triggerNames = event.triggerNames(*triggerResultsHandle).triggerNames();
  std::set<std::string> modules;
  for ( int i=0, n=triggerNames.size(); i<n; ++i ) {
    if ( !triggerResultsHandle->accept(i) ) continue;
    const auto& stmodules = hltConfig_.saveTagsModules(i);

    for ( size_t j=0, m=triggerPaths_.size(); j<m; ++j ) {
      const auto path = triggerPaths_[j];
      const auto module = (j < triggerModules_.size()) ? triggerModules_[j] : "";

      if ( hltConfig_.removeVersion(triggerNames[i]) != path ) continue;

      // Keep module names
      if ( module.empty() ) modules.insert(stmodules.begin(), stmodules.end()); // all modules if not specified
      else if ( std::find(stmodules.begin(), stmodules.end(), module) != stmodules.end() ) modules.insert(module);
    }
  }
  std::vector<math::XYZTLorentzVector> triggerObjectP4s;

  const auto& triggerObjects = triggerEventHandle->getObjects();
  for ( size_t keyIdx = 0; keyIdx < triggerEventHandle->sizeFilters(); ++keyIdx ) {
    if ( modules.count(triggerEventHandle->filterLabel(keyIdx)) == 0 ) continue;

    for ( auto objIdx : triggerEventHandle->filterKeys(keyIdx) ) {
      //if ( std::abs(triggerObjects[objIdx].id()) != 13 ) continue;
      triggerObjectP4s.push_back(triggerObjects[objIdx].particle().p4());
    }
  }
  if ( triggerObjectP4s.empty() ) return;

  // Select best tag muon
  for ( int i=0, n=muonHandle->size(); i<n; ++i ) {
    const auto& mu = muonHandle->at(i);
    const double pt = mu.pt();

    // Basic kinematic cuts
    if ( pt < minMuonPt_ or std::abs(mu.eta()) > maxMuonAbsEta_ ) continue;

    // Tight muon ID
    if ( mu.track().isNull() ) continue;
    if ( !muon::isTightMuon(mu, *goodPV) ) continue;

    // Tight muon isolation
    const double chIso = mu.pfIsolationR03().sumChargedHadronPt;
    const double nhIso = mu.pfIsolationR03().sumNeutralHadronEt;
    const double phIso = mu.pfIsolationR03().sumPhotonEt;
    const double puIso = mu.pfIsolationR03().sumPUPt;
    if ( chIso + std::max(0., nhIso+phIso-0.5*puIso) > pt*maxMuonRelIso_ ) continue;

    // Trigger matching
    const bool isTrigMatching = [&](){
      for ( const auto& to : triggerObjectP4s ) {
        if ( deltaR(mu, to) < 0.1 and std::abs(mu.pt()-to.pt()) < 0.5*mu.pt() ) return true;
      }
      return false; }();
    if ( !isTrigMatching ) return;

    if ( tagRef.isNull() or tagRef->pt() < pt ) tagRef = reco::MuonRef(muonHandle, i);
  }
  if ( tagRef.isNull() ) return;

  // Find best tag + probe pair
  for ( int i=0, n=muonHandle->size(); i<n; ++i ) {
    const auto& mu = muonHandle->at(i);
    if ( !mu.isTrackerMuon() ) continue;
    if ( !muon::isGoodMuon(mu, muon::TMOneStationLoose) ) continue;
    if ( mu.track()->originalAlgo() == reco::TrackBase::muonSeededStepOutIn ) continue;

    const double pt = mu.pt();

    // Basic kinematic cuts
    if ( pt < minTrackPt_ or std::abs(mu.eta()) > maxTrackAbsEta_ ) continue;

    // Require opposite charge and do overlap check
    if ( doCheckSign_ and mu.charge() == tagRef->charge() ) continue;
    if ( deltaR(mu, tagRef->p4()) < minDR_ ) continue;

    // Set four momentum with muon hypothesis, compute invariant mass
    const double m = (tagRef->p4()+mu.p4()).mass();
    if ( m < minMass_ or m > maxMass_ ) continue;

    if ( probeRef.isNull() or probeRef->pt() < pt ) {
      probeRef = reco::MuonRef(muonHandle, i);
      mass = m;
    }
  }
  if ( probeRef.isNull() ) return;

  // Now we have tag + probe pair, mass is already set to 'mass' variable
  h_tpMass_->Fill(mass);
  h_tagPt_->Fill(tagRef->pt());
  h_probePt_->Fill(probeRef->pt());

  // Now continue to find extrapolation from the probe track
  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  edm::Handle<RPCRecHitCollection> rpcHitHandle;
  event.getByToken(rpcHitToken_, rpcHitHandle);

  for ( auto match : probeRef->matches() ) {
    if ( match.detector() != MuonSubdetId::RPC ) continue;

    const RPCDetId rpcId(match.id);
    const RPCRoll* roll = rpcGeom->roll(match.id);
    const auto& bound = roll->surface().bounds();
    //const GlobalPoint gPos0 = roll->toGlobal(LocalPoint(0,0,0));

    LocalPoint lPos(match.x, match.y, 0);
    if ( !bound.inside(lPos) ) continue;

    const double match_x0 = roll->centreOfStrip(roll->strip(lPos)).x();

    //const auto gp = roll->toGlobal(lPos);
    const RPCDetId detId(match.id);

    bool isMatched = false;
    if ( detId.region() == 0 ) {
      const bool isInFiducial = (std::abs(lPos.y()) <= bound.length()/2-8 and
          std::abs(lPos.x()) <= bound.width()/2-8 );
      if ( !isInFiducial ) continue;

      // Find best-matching RPCRecHit
      auto rpcHitRange = rpcHitHandle->get(detId);
      if ( rpcHitRange.first != rpcHitRange.second ) {
        auto matchedHit = rpcHitRange.first;
        double minDX = 1e9;
        for ( auto rpcHit = rpcHitRange.first; rpcHit != rpcHitRange.second; ++rpcHit ) {
          const double dx = std::abs(rpcHit->localPosition().x()-match_x0);
          if ( dx < minDX ) {
            matchedHit = rpcHit;
            minDX = dx;
          }
        }
        //const auto hitLPos = matchedHit->localPosition();
        //const auto hitLErr = matchedHit->localPositionError();

        isMatched = true;
        /*const auto dp = hitLPos - lPos;
          vars[RESX] = dp.x();
          vars[RESY] = dp.y();
          vars[PULLX] = dp.x()/std::sqrt(hitLErr.xx());
          vars[PULLY] = dp.y()/std::sqrt(hitLErr.yy());
          vars[CLS] = matchedHit->clusterSize();
          vars[BX] = matchedHit->BunchX();*/
      }
    }
    else {
      const double wT = bound.width(), w0 = bound.widthAtHalfLength();
      const double slope = (wT-w0)/bound.length();
      const double w2AtY = slope*lPos.y() + w0/2;
      const bool isInFiducial = (std::abs(lPos.y()) <= bound.length()/2-8 and
          std::abs(lPos.x()) <= w2AtY-8 );
      if ( !isInFiducial ) continue;

      // Find best-matching RPCRecHit
      auto rpcHitRange = rpcHitHandle->get(detId);
      if ( rpcHitRange.first != rpcHitRange.second ) {
        auto matchedHit = rpcHitRange.first;
        double minDX = 1e9;
        for ( auto hit = rpcHitRange.first; hit != rpcHitRange.second; ++hit ) {
          //const double dr = std::hypot(hit->localPosition().x(), hit->localPosition().y());
          const double dx = std::abs(hit->localPosition().x()-match_x0);
          if ( dx < minDX ) {
            matchedHit = hit;
            minDX = dx;
          }
        }
        //const auto hitLPos = matchedHit->localPosition();
        //const auto hitLErr = matchedHit->localPositionError();

        isMatched = true;
        /*Iconst auto dp = hitLPos - lPos;
          vars[RESX] = dp.x();
          vars[RESY] = dp.y();
          vars[PULLX] = dp.x()/std::sqrt(hitLErr.xx());
          vars[PULLY] = dp.y()/std::sqrt(hitLErr.yy());
          vars[CLS] = matchedHit->clusterSize();
          vars[BX] = matchedHit->BunchX();*/
      }
    }

    const int bin = rawIdToBin_[match.id];
    h_detIdAll_->Fill(bin);
    if ( isMatched ) h_detIdPass_->Fill(bin);
  }
}

