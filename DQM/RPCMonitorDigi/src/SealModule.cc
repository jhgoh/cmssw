#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include <DQM/RPCMonitorDigi/interface/RPCFEDIntegrity.h>
DEFINE_FWK_MODULE(RPCFEDIntegrity);
#include <DQM/RPCMonitorDigi/interface/RPCMonitorRaw.h>
DEFINE_FWK_MODULE(RPCMonitorRaw);
#include <DQM/RPCMonitorDigi/interface/RPCMonitorLinkSynchro.h>
DEFINE_FWK_MODULE(RPCMonitorLinkSynchro);

#include "DQM/RPCMonitorDigi/interface/RPCMonitorDigi.h"
DEFINE_FWK_MODULE(RPCMonitorDigi);
#include "DQM/RPCMonitorDigi/interface/RPCRecHitProbability.h"
DEFINE_FWK_MODULE(RPCRecHitProbability);
#include "DQM/RPCMonitorDigi/interface/RPCTTUMonitor.h"
DEFINE_FWK_MODULE(RPCTTUMonitor);
