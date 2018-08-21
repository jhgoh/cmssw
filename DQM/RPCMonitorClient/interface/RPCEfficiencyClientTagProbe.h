#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

/*#include <DataFormats/MuonDetId/interface/RPCDetId.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include "FWCore/Framework/interface/ESHandle.h"

#include<string>
#include<map>
#include<fstream>
*/

class RPCEfficiencyClientTagProbe :public DQMEDHarvester
{
public:
  RPCEfficiencyClientTagProbe(const edm::ParameterSet&);
  ~RPCEfficiencyClientTagProbe() override = default;
  
protected:
  void dqmEndLuminosityBlock(DQMStore::IBooker&, DQMStore::IGetter&,
                             edm::LuminosityBlock const&, edm::EventSetup const&) override; //performed in the endLumi
  void dqmEndJob(DQMStore::IBooker&, DQMStore::IGetter&) override; //performed in the endJob

private:
  const unsigned int minNExtrapolation_;
  std::map<int, long> binToId_;

  typedef MonitorElement* MEP;

  MEP hEffBarrel_, hEffEndcap_;

};

