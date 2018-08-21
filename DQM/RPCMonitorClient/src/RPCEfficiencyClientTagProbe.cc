#include "DQM/RPCMonitorClient/interface/RPCEfficiencyClientTagProbe.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"

RPCEfficiencyClientTagProbe::RPCEfficiencyClientTagProbe(const edm::ParameterSet& pset):
  minNExtrapolation_(pset.getUntrackedParameter<unsigned int>("minNExtrapolation"))
{
}

void RPCEfficiencyClientTagProbe::dqmEndLuminosityBlock(DQMStore::IBooker& booker, DQMStore::IGetter& getter,
                                                        edm::LuminosityBlock const &, edm::EventSetup const& eventSetup)
{
  edm::ESHandle<RPCGeometry> rpcGeom;
  eventSetup.get<MuonGeometryRecord>().get(rpcGeom);
  auto hDetIdAll = getter.get("RPC/RPCEfficiency/TagAndProbe/detIdAll");
  auto axis = hDetIdAll->getTH1()->GetXaxis();

  for ( auto roll : rpcGeom->rolls() ) {
    const int rawId = roll->geographicalId().rawId();
    RPCDetId rpcId(rawId);
    const std::string name = RPCGeomServ(rpcId).name();
    const int bin = axis->FindBin(name.c_str());
    binToId_[bin] = rawId;
  }
}

void RPCEfficiencyClientTagProbe::dqmEndJob(DQMStore::IBooker& booker, DQMStore::IGetter& getter)
{
  booker.setCurrentFolder("RPC/RPCEfficiency");
  hEffBarrel_ = booker.book1D("EffBarrel", "Overall Efficiency Barrel;Efficiency (%);Number of rolls", 201, 0, 100.5);
  hEffEndcap_ = booker.book1D("EffEndcap", "Overall Efficiency Endcap;Efficiency (%);Number of rolls", 201, 0, 100.5);

  MEP hDetIdAll = getter.get("RPC/RPCEfficiency/TagAndProbe/detIdAll");
  MEP hDetIdPass = getter.get("RPC/RPCEfficiency/TagAndProbe/detIdPass");

  for ( int i=1, n=hDetIdAll->getNbinsX(); i<=n; ++i ) {
    const double countAll = hDetIdAll->getBinContent(i);
    const double countPass = hDetIdPass->getBinContent(i);

    if ( countAll < minNExtrapolation_ ) continue; // Skip to fill eff. with low statstics

    const int rawId = binToId_[i];
    RPCDetId detId(rawId);

    const double eff = countPass/countAll;
    if ( detId.region() == 0 ) hEffBarrel_->Fill(100*eff);
    else hEffEndcap_->Fill(100*eff);
  }

}
