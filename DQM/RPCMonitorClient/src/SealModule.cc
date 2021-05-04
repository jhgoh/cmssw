#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

//General Client
#include <DQM/RPCMonitorClient/interface/RPCDqmClient.h>
DEFINE_FWK_MODULE(RPCDqmClient);

#include <DQM/RPCMonitorClient/interface/RPCRecHitProbabilityClient.h>
DEFINE_FWK_MODULE(RPCRecHitProbabilityClient);

#include <DQM/RPCMonitorClient/interface/RPCChamberQuality.h>
DEFINE_FWK_MODULE(RPCChamberQuality);

#include <DQM/RPCMonitorClient/interface/RPCDcsInfoClient.h>
DEFINE_FWK_MODULE(RPCDcsInfoClient);

#include <DQM/RPCMonitorClient/interface/RPCEventSummary.h>
DEFINE_FWK_MODULE(RPCEventSummary);

