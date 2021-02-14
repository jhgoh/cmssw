#include "Validation/MuonRPCDigis/interface/RPCDigiPSetValid.h"

#include "DataFormats/Provenance/interface/StableProvenance.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include <TH1.h>

RPCDigiPSetValid::RPCDigiPSetValid(const edm::ParameterSet& pset):
  moduleLabel_(pset.getUntrackedParameter<std::string>("moduleLabel"))
{
}

void RPCDigiPSetValid::bookHistograms(DQMStore::IBooker& booker, edm::Run const& run, edm::EventSetup const& eventSetup)
{
  booker.setCurrentFolder("RPC/RPCDigisV/Info");

  //hNEvents_ = booker.book1D("NEvents", "Number of events;;Events", 6, 1, 7);
  // Area, nStrips, strip lengths, width...
  // Dead chambers, noisy chambers

  // SimHit entry angle (vs cluster size or number of digis from a muon)
  //hNSimHitPerStrip_ = booker.book1D("NSimHitPerStrip", "Number of simHits per strip;N simHits/strip;Entries", 100, 0, 

  isFilled_ = false;
  hDigiModel1_ = booker.bookString("DigiModel1", "");
  hDigiModel2_ = booker.bookString("DigiModel2", "");

  hPropSpeed1_ = booker.bookFloat("PropSpeed1");
  hPropSpeed2_ = booker.bookFloat("PropSpeed2");
  hTimeResol1_ = booker.bookFloat("TimeResol1");
  hTimeResol2_ = booker.bookFloat("TimeResol2");
  hTimeJitter1_ = booker.bookFloat("TimeJitter1");
  hTimeJitter2_ = booker.bookFloat("TimeJitter2");
  hClusterSize1_ = booker.bookFloat("ClusterSize1");
  hClusterSize2_ = booker.bookFloat("ClusterSize2");
  hEfficiency1_ = booker.bookFloat("Efficiency1");
  hEfficiency2_ = booker.bookFloat("Efficiency2");

  hIRPCTimeResol1_ = booker.bookFloat("IRPCTimeResol1");
  hIRPCTimeResol2_ = booker.bookFloat("IRPCTimeResol2");
  hIRPCTimeJitter1_ = booker.bookFloat("IRPCTimeJitter1");
  hIRPCTimeJitter2_ = booker.bookFloat("IRPCTimeJitter2");
}

void RPCDigiPSetValid::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
{
  using namespace std;

  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);

  //edm::Handle<edm::PSimHitContainer> simHitsHandle;
  //event.getByToken(simHitsToken_, simHitsHandle);

  if ( !isFilled_ ) {
    std::vector<std::string> procNames;

    std::vector<const edm::StableProvenance*> provenances;
    event.getAllStableProvenance(provenances);
    for ( const auto prov : provenances ) {
      const std::string label = prov->moduleLabel();
      const std::string instance = prov->productInstanceName();
      const std::string moduleLabel = instance.empty() ? label : label+":"+instance;

      if ( moduleLabel_ == moduleLabel ) {
        procNames.push_back(prov->processName());
      }
    }

    if ( !procNames.empty() ) {
      const std::string& procName = procNames.back();
      edm::ParameterSet pset;
      if ( event.getProcessParameterSet(procName, pset) ) {
        const auto modulePSet = pset.getParameter<edm::ParameterSet>(moduleLabel_);
        std::string tmpStr;
        tmpStr = modulePSet.getParameter<std::string>("digiModel");
        hDigiModel1_->Fill(tmpStr);
        tmpStr = modulePSet.getParameter<std::string>("digiIRPCModel");
        hDigiModel2_->Fill(tmpStr);
/*
  Noise: bool tracked  = true
  Signal: bool tracked  = true
  doBkgNoise: bool tracked  = false
*/

        const auto cfg1 = modulePSet.getParameter<edm::ParameterSet>("digiModelConfig");
        hPropSpeed1_->Fill(cfg1.getParameter<double>("signalPropagationSpeed"));
        hTimeResol1_->Fill(cfg1.getParameter<double>("timeResolution"));
        hTimeJitter1_->Fill(cfg1.getParameter<double>("timeJitter"));
        hEfficiency1_->Fill(cfg1.getParameter<double>("averageEfficiency"));
        hClusterSize1_->Fill(cfg1.getParameter<double>("averageClusterSize"));
/*
   BX_range: int32 tracked  = 5
   Frate: double tracked  = 1
   Gate: double tracked  = 25
   Nbxing: int32 tracked  = 9
   Rate: double tracked  = 0
   cosmics: bool tracked  = false
   deltatimeAdjacentStrip: double tracked  = 3
   digitizeElectrons: bool tracked  = true
   do_Y_coordinate: bool tracked  = true
   linkGateWidth: double tracked  = 20
   sigmaY: double tracked  = 2
   timingRPCOffset: double tracked  = 50
*/
        hIRPCTimeResol1_->Fill(cfg1.getParameter<double>("IRPC_time_resolution"));
        hIRPCTimeJitter1_->Fill(cfg1.getParameter<double>("IRPC_electronics_jitter"));

        const auto cfg2 = modulePSet.getParameter<edm::ParameterSet>("digiModelConfig");
        hPropSpeed2_->Fill(cfg2.getParameter<double>("signalPropagationSpeed"));
        hTimeResol2_->Fill(cfg2.getParameter<double>("timeResolution"));
        hTimeJitter2_->Fill(cfg2.getParameter<double>("timeJitter"));
        hEfficiency2_->Fill(cfg2.getParameter<double>("averageEfficiency"));
        hClusterSize2_->Fill(cfg2.getParameter<double>("averageClusterSize"));
        hIRPCTimeResol2_->Fill(cfg2.getParameter<double>("IRPC_time_resolution"));
        hIRPCTimeJitter2_->Fill(cfg2.getParameter<double>("IRPC_electronics_jitter"));
      }
    }

    isFilled_ = true;
  }
}
