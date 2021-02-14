#include "Validation/MuonRPCDigis/interface/RPCSimDigiValid.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include <TH1.h>

RPCSimDigiValid::RPCSimDigiValid(const edm::ParameterSet& pset)
{
  simHitsToken_ = consumes<edm::PSimHitContainer>(pset.getUntrackedParameter<edm::InputTag>("simHits"));
  rpcDigisToken_ = consumes<RPCDigiCollection>(pset.getUntrackedParameter<edm::InputTag>("rpcDigis"));
  genParticlesToken_ = consumes<reco::GenParticleCollection>(pset.getUntrackedParameter<edm::InputTag>("genParticles"));

}

void RPCSimDigiValid::bookHistograms(DQMStore::IBooker& booker, edm::Run const& run, edm::EventSetup const& eventSetup)
{
  booker.setCurrentFolder("RPC/RPCDigisV");

  hNEvents_ = booker.book1D("NEvents", "Number of events;;Events", 10, 1, 11);
  if ( TH1* h = hNEvents_->getTH1() ) {
    h->GetXaxis()->SetBinLabel(1, "All events");
    h->GetXaxis()->SetBinLabel(2, "has GenMuon");
    h->GetXaxis()->SetBinLabel(3, "GenMuon pt > 7 |#eta|<2.4");
    h->GetXaxis()->SetBinLabel(4, "GenMuon pt > 7 |#eta|<1.9");
    h->GetXaxis()->SetBinLabel(5, "has RPC simhit");
    h->GetXaxis()->SetBinLabel(6, "has RPC Digis");
  }

  // Area, nStrips, strip lengths, width...
  // Dead chambers, noisy chambers

  // Histograms for the simHits
  booker.setCurrentFolder("RPC/RPCDigisV/SimHits");

  hSimHitPid_ = booker.book1D("ParticleId", "SimHits Particle ID", 11, 1, 12);
  if ( TH1* h = hSimHitPid_->getTH1() ) {
    h->GetXaxis()->SetBinLabel(1, "#mu-");
    h->GetXaxis()->SetBinLabel(2, "#mu+");
    h->GetXaxis()->SetBinLabel(3, "e-");
    h->GetXaxis()->SetBinLabel(4, "e+");
    h->GetXaxis()->SetBinLabel(5, "#pi+");
    h->GetXaxis()->SetBinLabel(6, "#pi-");
    h->GetXaxis()->SetBinLabel(7, "K+");
    h->GetXaxis()->SetBinLabel(8, "K-");
    h->GetXaxis()->SetBinLabel(9, "p+");
    h->GetXaxis()->SetBinLabel(10, "p-");
    h->GetXaxis()->SetBinLabel(11, "Others");
  }

  // SimHit timing
  hSimHitTime_ = booker.book1D("SimHitTime", "SimHit time;Time (ns);Entries", 200, -30, 70);
  hSimHitTime_Mu_ = booker.book1D("SimHitTime_Mu", "SimHit time (#mu#pm);Time (ns);Entries", 200, -30, 70);
  hSimHitTime_El_ = booker.book1D("SimHitTime_El", "SimHit time (e#pm);Time (ns);Entries", 200, -30, 70);
  hSimHitTime_Others_ = booker.book1D("SimHitTime_Others", "SimHit time (others);Time (ns);Entries", 200, -30, 70);

  // SimHit energy loss
  const double maxELoss = 0.1;
  hSimHitELoss_ = booker.book1D("SimHitELoss", "SimHit energy loss;ELoss (keV);Entries", 100, 0, maxELoss);
  hSimHitELoss_Mu_ = booker.book1D("SimHitELoss_Mu", "SimHit energy loss (#mu#pm);ELoss (keV);Entries", 100, 0, maxELoss);
  hSimHitELoss_El_ = booker.book1D("SimHitELoss_El", "SimHit energy loss (e#pm);ELoss (keV);Entries", 100, 0, maxELoss);
  hSimHitELoss_Others_ = booker.book1D("SimHitELoss_Others", "SimHit energy loss (others);ELoss (keV);Entries", 100, 0, maxELoss);

  // Eta, phi, Z etc with respect to the particle IDs (muon electron or others)
  booker.setCurrentFolder("RPC/RPCDigisV/SimHits/Barrel");
  hSimHitELoss_Barrel_ = booker.book1D("SimHitELoss", "SimHit energy loss;ELoss (keV);Entries", 100, 0, maxELoss);
  hSimHitR_Barrel_ = booker.book1D("SimHitR", "SimHit distance from the beam line (Barrel);R (cm);Entries", 100, 400, 800);
  hSimHitZ_Barrel_ = booker.book1D("SimHitZ", "SimHit z position (Barrel);z (cm);Entries", 160, -800, 800);

  hSimHitZ_RB1in_  = booker.book1D("SimHitZ_RB1in" , "SimHit z position (RB1in);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB1out_ = booker.book1D("SimHitZ_RB1out", "SimHit z position (RB1out);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB2in_  = booker.book1D("SimHitZ_RB2in" , "SimHit z position (RB2in);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB2out_ = booker.book1D("SimHitZ_RB2out", "SimHit z position (RB2out);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB3_    = booker.book1D("SimHitZ_RB3"   , "SimHit z position (RB3);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB4_    = booker.book1D("SimHitZ_RB4"   , "SimHit z position (RB4);z (cm);Entries", 160, -800, 800);

  hSimHitZ_RB1in_Mu_  = booker.book1D("SimHitZ_RB1in_Mu" , "SimHit z position (RB1in, #mu#pm);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB1out_Mu_ = booker.book1D("SimHitZ_RB1out_Mu", "SimHit z position (RB1out, #mu#pm);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB2in_Mu_  = booker.book1D("SimHitZ_RB2in_Mu" , "SimHit z position (RB2in, #mu#pm);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB2out_Mu_ = booker.book1D("SimHitZ_RB2out_Mu", "SimHit z position (RB2out, #mu#pm);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB3_Mu_    = booker.book1D("SimHitZ_RB3_Mu"   , "SimHit z position (RB3, #mu#pm);z (cm);Entries", 160, -800, 800);
  hSimHitZ_RB4_Mu_    = booker.book1D("SimHitZ_RB4_Mu"   , "SimHit z position (RB4, #mu#pm);z (cm);Entries", 160, -800, 800);

  booker.setCurrentFolder("RPC/RPCDigisV/SimHits/Endcap+");
  hSimHitELoss_EndcapP_ = booker.book1D("SimHitELoss", "SimHit energy loss;ELoss (keV);Entries", 100, 0, maxELoss);
  hSimHitR_EndcapP_ = booker.book1D("SimHitR", "SimHit distance from the beam line (Endcap+);R (cm);Entries", 140, 100, 800);
  hSimHitZ_EndcapP_ = booker.book1D("SimHitAbsZ", "SimHit |z| (Endcap+);|z| (cm);Entries", 80, 650, 1050);

  hSimHitR_REP1_ = booker.book1D("SimHitR_RE+1", "SimHit distance from the beam line (RE+1);R (cm);Entries", 140, 100, 800);
  hSimHitR_REP2_ = booker.book1D("SimHitR_RE+2", "SimHit distance from the beam line (RE+2);R (cm);Entries", 140, 100, 800);
  hSimHitR_REP3_ = booker.book1D("SimHitR_RE+3", "SimHit distance from the beam line (RE+3);R (cm);Entries", 140, 100, 800);
  hSimHitR_REP4_ = booker.book1D("SimHitR_RE+4", "SimHit distance from the beam line (RE+4);R (cm);Entries", 140, 100, 800);

  hSimHitR_REP1_Mu_ = booker.book1D("SimHitR_RE+1_Mu", "SimHit distance from the beam line (RE+1, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REP2_Mu_ = booker.book1D("SimHitR_RE+2_Mu", "SimHit distance from the beam line (RE+2, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REP3_Mu_ = booker.book1D("SimHitR_RE+3_Mu", "SimHit distance from the beam line (RE+3, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REP4_Mu_ = booker.book1D("SimHitR_RE+4_Mu", "SimHit distance from the beam line (RE+4, #mu#pm);R (cm);Entries", 140, 100, 800);

  booker.setCurrentFolder("RPC/RPCDigisV/SimHits/Endcap-");
  hSimHitELoss_EndcapM_ = booker.book1D("SimHitELoss", "SimHit energy loss;ELoss (keV);Entries", 100, 0, maxELoss);
  hSimHitR_EndcapM_ = booker.book1D("SimHitR", "SimHit distance from the beam line (Endcap-);R (cm);Entries", 140, 100, 800);
  hSimHitZ_EndcapM_ = booker.book1D("SimHitAbsZ", "SimHit |z| (Endcap-);|z| (cm);Entries", 80, 650, 1050);

  hSimHitR_REM1_ = booker.book1D("SimHitR_RE-1", "SimHit distance from the beam line (RE-1, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REM2_ = booker.book1D("SimHitR_RE-2", "SimHit distance from the beam line (RE-2, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REM3_ = booker.book1D("SimHitR_RE-3", "SimHit distance from the beam line (RE-3, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REM4_ = booker.book1D("SimHitR_RE-4", "SimHit distance from the beam line (RE-4, #mu#pm);R (cm);Entries", 140, 100, 800);

  hSimHitR_REM1_Mu_ = booker.book1D("SimHitR_RE-1_Mu", "SimHit distance from the beam line (RE-1, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REM2_Mu_ = booker.book1D("SimHitR_RE-2_Mu", "SimHit distance from the beam line (RE-2, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REM3_Mu_ = booker.book1D("SimHitR_RE-3_Mu", "SimHit distance from the beam line (RE-3, #mu#pm);R (cm);Entries", 140, 100, 800);
  hSimHitR_REM4_Mu_ = booker.book1D("SimHitR_RE-4_Mu", "SimHit distance from the beam line (RE-4, #mu#pm);R (cm);Entries", 140, 100, 800);

  // SimHit entry angle (vs cluster size or number of digis from a muon)
  //hNSimHitPerStrip_ = booker.book1D("NSimHitPerStrip", "Number of simHits per strip;N simHits/strip;Entries", 100, 0, 
}

void RPCSimDigiValid::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  edm::Handle<edm::PSimHitContainer> simHitsHandle;
  event.getByToken(simHitsToken_, simHitsHandle);

  edm::Handle<RPCDigiCollection> rpcDigisHandle;
  event.getByToken(rpcDigisToken_, rpcDigisHandle);

  edm::Handle<reco::GenParticleCollection> genParticlesHandle;
  event.getByToken(genParticlesToken_, genParticlesHandle);

  hNEvents_->Fill(1);

  // Check generator level muons in RPC acceptance
  int nGenMuon = 0, nGenMuonEta24 = 0, nGenMuonEta19 = 0;
  if ( genParticlesHandle.isValid() ) {
    for ( auto p : *genParticlesHandle ) {
      if ( p.status() != 1 or std::abs(p.pdgId()) != 13 ) continue;
      ++nGenMuon;
      if ( p.pt() < 7 ) continue;
      const double aeta = std::abs(p.eta());
      if ( aeta < 2.5 ) ++nGenMuonEta24;
      if ( aeta < 1.9 ) ++nGenMuonEta19;
    }
  }
  if ( nGenMuon > 0 ) hNEvents_->Fill(2);
  if ( nGenMuonEta24 > 0 ) hNEvents_->Fill(3);
  if ( nGenMuonEta19 > 0 ) hNEvents_->Fill(4);

  // Fill g4SimHits plots
  int nRPCSimHits = 0;
  if ( simHitsHandle.isValid() ) {
    for ( auto& simHit : *simHitsHandle ) {
      const RPCDetId detId(simHit.detUnitId());
      const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(detId));
      if ( !roll ) continue;
      ++nRPCSimHits;

      // Particle types
      const int pid = simHit.particleType();
      switch ( pid ) {
        case   +13: hSimHitPid_->Fill( 1); break;
        case   -13: hSimHitPid_->Fill( 2); break;
        case   +11: hSimHitPid_->Fill( 3); break;
        case   -11: hSimHitPid_->Fill( 4); break;
        case  +211: hSimHitPid_->Fill( 5); break;
        case  -211: hSimHitPid_->Fill( 6); break;
        case  +321: hSimHitPid_->Fill( 7); break;
        case  -321: hSimHitPid_->Fill( 8); break;
        case +2212: hSimHitPid_->Fill( 9); break;
        case -2212: hSimHitPid_->Fill(10); break;
        default:    hSimHitPid_->Fill(11);
      }

      // SimHit position in the global coordinate
      const auto simGPos = roll->toGlobal(simHit.localPosition());
      const double r = simGPos.perp();
      const double z = simGPos.z();

      // SimHit ELoss and Time
      const auto refGPos = roll->toGlobal(LocalPoint(simHit.localPosition().x(),0,0));
      const double time = simHit.timeOfFlight() - refGPos.mag()/speedOfLight_;
      const double eLoss = simHit.energyLoss()*1e3; // MeV to keV

      if ( detId.region() == 0 ) {
        hSimHitELoss_Barrel_->Fill(eLoss);
        hSimHitR_Barrel_->Fill(r);
        hSimHitZ_Barrel_->Fill(z);

        const int station = detId.station();
        const int layer = detId.layer();
        if ( station == 1 and layer == 1 ) hSimHitZ_RB1in_->Fill(z);
        else if ( station == 1 and layer == 2 ) hSimHitZ_RB1out_->Fill(z);
        else if ( station == 2 and layer == 1 ) hSimHitZ_RB2in_->Fill(z);
        else if ( station == 2 and layer == 2 ) hSimHitZ_RB2out_->Fill(z);
        else if ( station == 3 ) hSimHitZ_RB3_->Fill(z);
        else if ( station == 4 ) hSimHitZ_RB4_->Fill(z);

        if ( std::abs(pid) == 13 ) {
          if ( station == 1 and layer == 1 ) hSimHitZ_RB1in_Mu_->Fill(z);
          else if ( station == 1 and layer == 2 ) hSimHitZ_RB1out_Mu_->Fill(z);
          else if ( station == 2 and layer == 1 ) hSimHitZ_RB2in_Mu_->Fill(z);
          else if ( station == 2 and layer == 2 ) hSimHitZ_RB2out_Mu_->Fill(z);
          else if ( station == 3 ) hSimHitZ_RB3_Mu_->Fill(z);
          else if ( station == 4 ) hSimHitZ_RB4_Mu_->Fill(z);
        }
      }
      else if ( detId.region() == 1 ) {
        hSimHitELoss_EndcapP_->Fill(eLoss);
        hSimHitR_EndcapP_->Fill(r);
        hSimHitZ_EndcapP_->Fill(z);

        const int station = detId.station();
        if ( station == 1 ) hSimHitR_REP1_->Fill(r);
        else if ( station == 2 ) hSimHitR_REP2_->Fill(r);
        else if ( station == 3 ) hSimHitR_REP3_->Fill(r);
        else if ( station == 4 ) hSimHitR_REP4_->Fill(r);

        if ( std::abs(pid) == 13 ) {
          if ( station == 1 ) hSimHitR_REP1_Mu_->Fill(r);
          else if ( station == 2 ) hSimHitR_REP2_Mu_->Fill(r);
          else if ( station == 3 ) hSimHitR_REP3_Mu_->Fill(r);
          else if ( station == 4 ) hSimHitR_REP4_Mu_->Fill(r);
        }
      }
      else if ( detId.region() == -1 ) {
        hSimHitELoss_EndcapM_->Fill(eLoss);
        hSimHitR_EndcapM_->Fill(r);
        hSimHitZ_EndcapM_->Fill(-z);

        const int station = detId.station();
        if ( station == 1 ) hSimHitR_REM1_->Fill(r);
        else if ( station == 2 ) hSimHitR_REM2_->Fill(r);
        else if ( station == 3 ) hSimHitR_REM3_->Fill(r);
        else if ( station == 4 ) hSimHitR_REM4_->Fill(r);

        if ( std::abs(pid) == 13 ) {
          if ( station == 1 ) hSimHitR_REM1_Mu_->Fill(r);
          else if ( station == 2 ) hSimHitR_REM2_Mu_->Fill(r);
          else if ( station == 3 ) hSimHitR_REM3_Mu_->Fill(r);
          else if ( station == 4 ) hSimHitR_REM4_Mu_->Fill(r);
        }
      }

      hSimHitTime_->Fill(time);
      hSimHitELoss_->Fill(eLoss);

      if ( std::abs(pid) == 13 ) {
        hSimHitTime_Mu_->Fill(time);
        hSimHitELoss_Mu_->Fill(eLoss);
      }
      else if ( std::abs(pid) == 11 ) {
        hSimHitTime_El_->Fill(time);
        hSimHitELoss_El_->Fill(eLoss);
      }
      else {
        hSimHitTime_Others_->Fill(time);
        hSimHitELoss_Others_->Fill(eLoss);
      }

    }
  }
  if ( nRPCSimHits > 0 ) hNEvents_->Fill(5);

  // Fill rpcDigis plots
  int nRPCDigis = 0;
  for ( auto rpcDigisItr = rpcDigisHandle->begin(); rpcDigisItr != rpcDigisHandle->end(); ++rpcDigisItr ) {
    const RPCDetId detId((*rpcDigisItr).first);
    const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeom->roll(detId));
    if ( !roll ) continue;
    const RPCDigiCollection::Range& range = (*rpcDigisItr).second;
    const int nDigis = range.second-range.first+1;
    nRPCDigis += nDigis;

    //for ( auto digiItr = range.first; digiItr != range.second; ++digiItr ) {
    //  
    //}
  }
  if ( nRPCDigis > 0 ) hNEvents_->Fill(6);
}
