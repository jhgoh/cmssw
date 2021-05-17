#include "FWCore/Framework/interface/MakerMacros.h"

//General Client
#include <DQM/RPCMonitorClient/plugins/RPCDqmClient.h>
DEFINE_FWK_MODULE(RPCDqmClient);

#include <DQM/RPCMonitorClient/plugins/RPCRecHitProbabilityClient.h>
DEFINE_FWK_MODULE(RPCRecHitProbabilityClient);

#include <DQM/RPCMonitorClient/plugins/RPCChamberQuality.h>
DEFINE_FWK_MODULE(RPCChamberQuality);

#include <DQM/RPCMonitorClient/plugins/RPCDcsInfoClient.h>
DEFINE_FWK_MODULE(RPCDcsInfoClient);

#include <DQM/RPCMonitorClient/plugins/RPCEventSummary.h>
DEFINE_FWK_MODULE(RPCEventSummary);

#include <DQM/RPCMonitorClient/plugins/RPCDaqInfo.h>
DEFINE_FWK_MODULE(RPCDaqInfo);

#include <DQM/RPCMonitorClient/plugins/RPCDCSSummary.h>
DEFINE_FWK_MODULE(RPCDCSSummary);

#include <DQM/RPCMonitorClient/plugins/RPCDataCertification.h>
DEFINE_FWK_MODULE(RPCDataCertification);
