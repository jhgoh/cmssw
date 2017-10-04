#ifndef EventFilter_RPCRawToDigi_RPCTwinMuxUnpacker_h
#define EventFilter_RPCRawToDigi_RPCTwinMuxUnpacker_h

#include <map>
#include <utility>
#include <set>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "CondFormats/DataRecord/interface/RPCTwinMuxLinkMapRcd.h"
#include "CondFormats/RPCObjects/interface/RPCAMCLinkMap.h"
#include "CondFormats/RPCObjects/interface/RPCLBLinkMap.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "EventFilter/RPCRawToDigi/interface/RPCAMC13Record.h"
#include "EventFilter/RPCRawToDigi/interface/RPCTwinMuxRecord.h"
#include "EventFilter/RPCRawToDigi/plugins/RPCAMCUnpacker.h"

class RPCAMCLinkCounters;

class RPCTwinMuxUnpacker
    : public RPCAMCUnpacker
{
public:
    RPCTwinMuxUnpacker(edm::stream::EDProducerBase & producer, edm::ParameterSet const & config);

    void beginRun(edm::Run const & run, edm::EventSetup const & setup);
    void produce(edm::Event & event, edm::EventSetup const & setup
                 , std::map<RPCAMCLink, rpcamc13::AMCPayload> const & amc_payload);

protected:
    bool processTwinMux(RPCAMCLink const & link, rpcamc13::AMCPayload const & payload
                        , RPCAMCLinkCounters & counters
                        , std::set<std::pair<RPCDetId, RPCDigi> > & digis) const;
    void processRPCRecord(RPCAMCLink link
                          , unsigned int bx_counter
                          , rpctwinmux::RPCRecord const & record
                          , RPCAMCLinkCounters & counters
                          , std::set<std::pair<RPCDetId, RPCDigi> > & digis
                          , int bx_min, int bx_max
                          , unsigned int input, unsigned int input_max) const;
    void putRPCDigis(edm::Event & event
                     , std::set<std::pair<RPCDetId, RPCDigi> > const & digis) const;

protected:
    bool fill_counters_;
    int bx_min_, bx_max_;

    edm::ESWatcher<RPCTwinMuxLinkMapRcd> es_tm_link_map_watcher_;
    edm::ESHandle<RPCAMCLinkMap> es_tm_link_map_;
    edm::ESHandle<RPCLBLinkMap> es_lb_link_map_;
};

#endif // EventFilter_RPCRawToDigi_RPCTwinMuxUnpacker_h
