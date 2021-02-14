#include "Validation/MuonRPCDigis/interface/RPCDigiConfigValid.h"

#include "DataFormats/Provenance/interface/StableProvenance.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include <TH1.h>

RPCDigiConfigValid::RPCDigiConfigValid(const edm::ParameterSet& pset):
  moduleLabel_(pset.getUntrackedParameter<std::string>("moduleLabel"))
{
}

void RPCDigiConfigValid::bookHistograms(DQMStore::IBooker& booker, edm::Run const& run, edm::EventSetup const& eventSetup)
{
  booker.setCurrentFolder("RPC/RPCDigisV/Config");

  isFilled_ = false;
  hDigiModel1_ = booker.bookString("DigiModel1", "");
  hDigiModel2_ = booker.bookString("DigiModel2", "");
  hNoise_ = booker.bookInt("Noise");
  hSignal_ = booker.bookInt("Signal");
  hDoBkgNoise_ = booker.bookInt("DoBkgNoise");

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

  hDigiElectron1_ = booker.bookInt("DigiElectron1");
  hDigiElectron2_ = booker.bookInt("DigiElectron2");
  hLinkGateWidth1_ = booker.bookFloat("LinkGateWidth1");
  hLinkGateWidth2_ = booker.bookFloat("LinkGateWidth2");

  hIRPCTimeResol1_ = booker.bookFloat("IRPCTimeResol1");
  hIRPCTimeResol2_ = booker.bookFloat("IRPCTimeResol2");
  hIRPCTimeJitter1_ = booker.bookFloat("IRPCTimeJitter1");
  hIRPCTimeJitter2_ = booker.bookFloat("IRPCTimeJitter2");

  // Area, nStrips, strip lengths, width...
  // Dead chambers, noisy chambers
}

void RPCDigiConfigValid::analyze(const edm::Event& event, const edm::EventSetup& eventSetup)
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
      edm::ParameterSet cfg;
      if ( event.getProcessParameterSet(procName, cfg) ) {
        const auto pset = cfg.getParameter<edm::ParameterSet>(moduleLabel_);
        std::string tmpStr;
        tmpStr = pset.getParameter<std::string>("digiModel");
        hDigiModel1_->Fill(tmpStr);
        tmpStr = pset.getParameter<std::string>("digiIRPCModel");
        hDigiModel2_->Fill(tmpStr);
        hNoise_->Fill(pset.getParameter<bool>("Noise"));
        hSignal_->Fill(pset.getParameter<bool>("Signal"));
        hDoBkgNoise_->Fill(pset.getParameter<bool>("doBkgNoise"));

        const auto pset1 = pset.getParameter<edm::ParameterSet>("digiModelConfig");
        hPropSpeed1_->Fill(pset1.getParameter<double>("signalPropagationSpeed"));
        hTimeResol1_->Fill(pset1.getParameter<double>("timeResolution"));
        hTimeJitter1_->Fill(pset1.getParameter<double>("timeJitter"));
        hEfficiency1_->Fill(pset1.getParameter<double>("averageEfficiency"));
        hClusterSize1_->Fill(pset1.getParameter<double>("averageClusterSize"));
        hDigiElectron1_->Fill(pset1.getParameter<bool>("digitizeElectrons"));
        hLinkGateWidth1_->Fill(pset1.getParameter<double>("linkGateWidth"));
        hIRPCTimeResol1_->Fill(pset1.getParameter<double>("IRPC_time_resolution"));
        hIRPCTimeJitter1_->Fill(pset1.getParameter<double>("IRPC_electronics_jitter"));

        const auto pset2 = pset.getParameter<edm::ParameterSet>("digiModelConfig");
        hPropSpeed2_->Fill(pset2.getParameter<double>("signalPropagationSpeed"));
        hTimeResol2_->Fill(pset2.getParameter<double>("timeResolution"));
        hTimeJitter2_->Fill(pset2.getParameter<double>("timeJitter"));
        hEfficiency2_->Fill(pset2.getParameter<double>("averageEfficiency"));
        hClusterSize2_->Fill(pset2.getParameter<double>("averageClusterSize"));
        hDigiElectron2_->Fill(pset2.getParameter<bool>("digitizeElectrons"));
        hLinkGateWidth2_->Fill(pset2.getParameter<double>("linkGateWidth"));
        hIRPCTimeResol2_->Fill(pset2.getParameter<double>("IRPC_time_resolution"));
        hIRPCTimeJitter2_->Fill(pset2.getParameter<double>("IRPC_electronics_jitter"));
      }
    }

    isFilled_ = true;
  }
}
