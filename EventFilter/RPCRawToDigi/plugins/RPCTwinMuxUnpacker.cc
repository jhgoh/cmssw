#include "EventFilter/RPCRawToDigi/plugins/RPCTwinMuxUnpacker.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/stream/EDProducerBase.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"

#include "CondFormats/DataRecord/interface/RPCLBLinkMapRcd.h"
#include "CondFormats/RPCObjects/interface/RPCAMCLink.h"
#include "DataFormats/RPCDigi/interface/RPCAMCLinkCounters.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "EventFilter/RPCRawToDigi/interface/RPCAMC13Record.h"
#include "EventFilter/RPCRawToDigi/interface/RPCAMCLinkEvents.h"
#include "EventFilter/RPCRawToDigi/interface/RPCAMCRecord.h"
#include "EventFilter/RPCRawToDigi/interface/RPCTwinMuxRecord.h"

RPCTwinMuxUnpacker::RPCTwinMuxUnpacker(edm::stream::EDProducerBase & producer, edm::ParameterSet const & config)
    : RPCAMCUnpacker(producer, config)
    , fill_counters_(config.getParameter<bool>("fillAMCCounters"))
    , bx_min_(config.getParameter<int>("bxMin"))
    , bx_max_(config.getParameter<int>("bxMax"))
{
    producer.produces<RPCDigiCollection>();
    if (fill_counters_) {
        producer.produces<RPCAMCLinkCounters>("RPCAMCUnpacker");
    }
}

void RPCTwinMuxUnpacker::beginRun(edm::Run const & run, edm::EventSetup const & setup)
{
    if (es_tm_link_map_watcher_.check(setup)) {
        setup.get<RPCTwinMuxLinkMapRcd>().get(es_tm_link_map_);
        std::set<int> feds;
        for (auto const & tm_link : es_tm_link_map_->getMap()) {
            feds.insert(tm_link.first.getFED());
        }
        feds_.assign(feds.begin(), feds.end());
    }
}

void RPCTwinMuxUnpacker::produce(edm::Event & event, edm::EventSetup const & setup
                                 , std::map<RPCAMCLink, rpcamc13::AMCPayload> const & amc_payload)
{
    // Get EventSetup Electronics Maps
    setup.get<RPCTwinMuxLinkMapRcd>().get(es_tm_link_map_);
    setup.get<RPCLBLinkMapRcd>().get(es_lb_link_map_);

    std::set<std::pair<RPCDetId, RPCDigi> > rpc_digis;
    std::unique_ptr<RPCAMCLinkCounters> counters(new RPCAMCLinkCounters());

    for (std::pair<RPCAMCLink const, rpcamc13::AMCPayload> const & payload : amc_payload) {
        processTwinMux(payload.first, payload.second
                       , *counters, rpc_digis);
    }

    putRPCDigis(event, rpc_digis);
    if (fill_counters_) {
        event.put(std::move(counters), "RPCAMCUnpacker");
    }
}

bool RPCTwinMuxUnpacker::processTwinMux(RPCAMCLink const & link, rpcamc13::AMCPayload const & payload
                                        , RPCAMCLinkCounters & counters
                                        , std::set<std::pair<RPCDetId, RPCDigi> > & digis) const
{
    LogDebug("RPCTwinMuxRawToDigi") << "TwinMux " << link
                                    << ", size " << payload.getData().size();

    if (!payload.isValid()) {
        return false;
    }

    std::vector<std::uint64_t>::const_iterator word(payload.getData().begin());
    rpctwinmux::Header header(&(*word));
    ++word; ++word;
    std::vector<std::uint64_t>::const_iterator word_end(payload.getData().end());
    --word_end; // rpcamc::Trailer
    unsigned int bx_counter(header.getBXCounter());

    int bx_min(bx_min_), bx_max(bx_max_);
    if (header.hasRPCBXWindow()) {
        bx_min = std::max(bx_min, header.getRPCBXMin());
        bx_max = std::min(bx_max, header.getRPCBXMax());
        LogDebug("RPCTwinMuxRawToDigi") << "BX range set to " << bx_min << ", " << bx_max;
    }

    bool has_first_rpc_word(false);
    rpctwinmux::RPCRecord rpc_record;
    for ( ; word != word_end ; ++word) {
        unsigned int type(rpctwinmux::Record::getType(*word));
        LogDebug("RPCTwinMuxRawToDigi") << "TwinMux data type " << std::hex << type << std::dec;
        if (type == rpctwinmux::Record::rpc_first_type_) {
            if (has_first_rpc_word) {
                processRPCRecord(link, bx_counter, rpc_record, counters, digis, bx_min, bx_max, 0, 1);
            }
            rpc_record.reset();
            rpc_record.set(0, *word);
            has_first_rpc_word = true;
        } else if (type == rpctwinmux::Record::rpc_second_type_) {
            if (!has_first_rpc_word) {
                edm::LogWarning("RPCTwinMuxRawToDigi") << "Received second RPC word without first";
            } else {
                rpc_record.set(1, *word);
                processRPCRecord(link, bx_counter, rpc_record, counters, digis, bx_min, bx_max, 0, 4);
                has_first_rpc_word = false;
            }
        }
    }
    if (has_first_rpc_word) {
        processRPCRecord(link, bx_counter, rpc_record, counters, digis, bx_min, bx_max, 0, 1);
    }

    return true;
}

void RPCTwinMuxUnpacker::processRPCRecord(RPCAMCLink link
                                          , unsigned int bx_counter
                                          , rpctwinmux::RPCRecord const & record
                                          , RPCAMCLinkCounters & counters
                                          , std::set<std::pair<RPCDetId, RPCDigi> > & digis
                                          , int bx_min, int bx_max
                                          , unsigned int input, unsigned int input_max) const
{
    LogDebug("RPCTwinMuxRawToDigi") << "RPCRecord " << std::hex << record.getRecord()[0]
                                    << ", " << record.getRecord()[1] << std::dec << std::endl;
    unsigned int fed(link.getFED());
    unsigned int amc_number(link.getAMCNumber());
    int bx_offset(record.getBXOffset());
    for ( ; input <= input_max ; ++input) {
        link.setAMCInput(input);
        rpctwinmux::RPCLinkRecord link_record(record.getRPCLinkRecord(input));

        if (link_record.isError()) {
            if (fill_counters_ && bx_offset == 0) {
                counters.add(RPCAMCLinkEvents::input_link_error_, link);
            }
            LogDebug("RPCTwinMuxRawToDigi") << "Link in error for " << link;
            continue;
        } else if (!link_record.isAcknowledge()) {
            if (fill_counters_ && bx_offset == 0) {
                counters.add(RPCAMCLinkEvents::input_link_ack_fail_, link);
            }
            LogDebug("RPCTwinMuxRawToDigi") << "Link without acknowledge for " << link;
            continue;
        }

        std::uint8_t data(link_record.getPartitionData());
        if (!data) {
            continue;
        }

        int bx(bx_offset - (int)(link_record.getDelay()));
        LogDebug("RPCTwinMuxRawToDigi") << "RPC BX " << bx << " for offset " << bx_offset;

        if (fill_counters_ && bx == 0 && link_record.isEOD()) { // EOD comes at the last delay
            counters.add(RPCAMCLinkEvents::input_eod_, link);
        }

        RPCAMCLinkMap::map_type::const_iterator link_it = es_tm_link_map_->getMap().find(link);
        if (link_it == es_tm_link_map_->getMap().end()) {
            if (fill_counters_ && bx_offset == 0) {
                counters.add(RPCAMCLinkEvents::amc_link_invalid_, RPCAMCLink(fed, amc_number));
            }
            LogDebug("RPCTwinMuxRawToDigi") << "Skipping unknown TwinMuxLink " << link;
            continue;
        }

        RPCLBLink lb_link = link_it->second;

        if (link_record.getLinkBoard() > (unsigned int)RPCLBLink::max_linkboard_) {
            if (fill_counters_ && bx_offset == 0) {
                counters.add(RPCAMCLinkEvents::input_lb_invalid_, link);
            }
            LogDebug("RPCTwinMuxRawToDigi") << "Skipping invalid LinkBoard " << link_record.getLinkBoard()
                                            << " for record " << link << " (" << std::hex << link_record.getRecord()
                                            << " in " << record.getRecord()[0] << ':' << record.getRecord()[1] << std::dec
                                            << " from " << link;
            continue;
        }

        if (link_record.getConnector() > (unsigned int)RPCLBLink::max_connector_) {
            if (fill_counters_ && bx_offset == 0) {
                counters.add(RPCAMCLinkEvents::input_connector_invalid_, link);
            }
            LogDebug("RPCTwinMuxRawToDigi") << "Skipping invalid Connector " << link_record.getConnector()
                                            << " for record " << link << " (" << std::hex << link_record.getRecord()
                                            << " in " << record.getRecord()[0] << ':' << record.getRecord()[1] << std::dec
                                            << ") from " << link;
            continue;
        }

        lb_link.setLinkBoard(link_record.getLinkBoard());
        lb_link.setConnector(link_record.getConnector());

        RPCLBLinkMap::map_type::const_iterator lb_link_it = es_lb_link_map_->getMap().find(lb_link);
        if (lb_link_it == es_lb_link_map_->getMap().end()) {
            if (fill_counters_ && bx_offset == 0) {
                counters.add(RPCAMCLinkEvents::input_connector_not_used_, link);
            }
            LogDebug("RPCTwinMuxRawToDigi") << "Could not find " << lb_link
                                            << " for record " << link << " (" << std::hex << link_record.getRecord()
                                            << " in " << record.getRecord()[0] << ':' << record.getRecord()[1] << std::dec
                                            << ") from " << link;
            continue;
        }

        if (bx < bx_min || bx > bx_max) {
            continue;
        }

        if (fill_counters_ && bx == 0) {
            counters.add(RPCAMCLinkEvents::input_event_, link);
        }

        RPCFebConnector const & feb_connector(lb_link_it->second);
        RPCDetId det_id(feb_connector.getRPCDetId());
        unsigned int channel_offset(link_record.getPartition() ? 9 : 1); // 1-16

        for (unsigned int channel = 0 ; channel < 8 ; ++channel) {
            if (data & (0x1 << channel)) {
                unsigned int strip(feb_connector.getStrip(channel + channel_offset));
                if (strip) {
                    digis.insert(std::pair<RPCDetId, RPCDigi>(det_id, RPCDigi(strip, bx)));
                    LogDebug("RPCTwinMuxRawToDigi") << "RPCDigi " << det_id.rawId()
                                                    << ", " << strip << ", " << bx;
                }
            }
        }

        // rpctwinmux::RPCBXRecord checks postponed: not implemented in firmware as planned and tbd if design or firmware should change
    }
}

void RPCTwinMuxUnpacker::putRPCDigis(edm::Event & event
                                     , std::set<std::pair<RPCDetId, RPCDigi> > const & digis) const
{
    std::unique_ptr<RPCDigiCollection> rpc_digi_collection(new RPCDigiCollection());
    RPCDetId rpc_det_id;
    std::vector<RPCDigi> local_digis;
    for (std::pair<RPCDetId, RPCDigi> const & rpc_digi : digis) {
        LogDebug("RPCTwinMuxRawToDigi") << "RPCDigi " << rpc_digi.first.rawId()
                                        << ", " << rpc_digi.second.strip() << ", " << rpc_digi.second.bx();
        if (rpc_digi.first != rpc_det_id) {
            if (!local_digis.empty()) {
                rpc_digi_collection->put(RPCDigiCollection::Range(local_digis.begin(), local_digis.end()), rpc_det_id);
                local_digis.clear();
            }
            rpc_det_id = rpc_digi.first;
        }
        local_digis.push_back(rpc_digi.second);
    }
    if (!local_digis.empty()) {
        rpc_digi_collection->put(RPCDigiCollection::Range(local_digis.begin(), local_digis.end()), rpc_det_id);
    }

    event.put(std::move(rpc_digi_collection));
}


#include "FWCore/Framework/interface/MakerMacros.h"
#include "EventFilter/RPCRawToDigi/plugins/RPCAMCUnpackerFactory.h"
DEFINE_EDM_PLUGIN (RPCAMCUnpackerFactory, RPCTwinMuxUnpacker, "RPCTwinMuxUnpacker");
