#include "EventFilter/RPCRawToDigi/interface/RPCTwinMuxRecord.h"

namespace rpctwinmux {

Header::Header()
{
    record_[0] = 0x0;
    record_[1] = dt_bx_window_mask_ | rpc_bx_window_mask_ | ho_bx_window_mask_;
}

Header::Header(std::uint64_t const record[2])
    : rpcamc::Header(record)
{}

Header::Header(unsigned int amc_number
               , unsigned int event_counter
               , unsigned int bx_counter
               , unsigned int data_length
               , unsigned int orbit_counter
               , unsigned int board_id
               , unsigned int dt_bx_window
               , int rpc_bx_min, int rpc_bx_max
               , unsigned int ho_bx_window)
    : rpcamc::Header(amc_number, event_counter, bx_counter, data_length, orbit_counter, board_id)
{
    setDTBXWindow(dt_bx_window);
    setRPCBXWindow(rpc_bx_min, rpc_bx_max);
    setHOBXWindow(ho_bx_window);
}

Header::~Header()
{}

Record::Record(std::uint64_t const record)
    : record_(record)
{}

RPCLinkRecord::RPCLinkRecord(std::uint32_t const record)
    : record_(record)
{}

RPCBXRecord::RPCBXRecord(std::uint8_t const record)
    : record_(record)
{}

unsigned int const RPCRecord::link_record_word_[] = {0, 0, 1, 1, 1};
unsigned int const RPCRecord::link_record_offset_[] = {20, 0, 40, 20, 0};
unsigned int const RPCRecord::bx_record_offset_[] = {52, 49, 46, 43, 40};

RPCRecord::RPCRecord()
{
    record_[0] = Record::rpc_first_identifier_
        | ((std::uint64_t)RPCLinkRecord::da_mask_ << link_record_offset_[0])
        | ((std::uint64_t)RPCLinkRecord::da_mask_ << link_record_offset_[1]);
    record_[1] = Record::rpc_second_identifier_
        | ((std::uint64_t)RPCLinkRecord::da_mask_ << link_record_offset_[2])
        | ((std::uint64_t)RPCLinkRecord::da_mask_ << link_record_offset_[3])
        | ((std::uint64_t)RPCLinkRecord::da_mask_ << link_record_offset_[4]);
}

RPCRecord::RPCRecord(std::uint64_t const record[2])
{
    record_[0] = record[0];
    record_[1] = record[1];
}

} // namespace rpctwinmux
