#ifndef EventFilter_RPCRawToDigi_RPCTwinMuxRecord_h
#define EventFilter_RPCRawToDigi_RPCTwinMuxRecord_h

#include <cstdint>

#include "EventFilter/RPCRawToDigi/interface/RPCAMCRecord.h"

namespace rpctwinmux {

class Header
    : public rpcamc::Header
{
protected:
    // Second word
    static std::uint64_t const dt_bx_window_mask_  = 0x0000001f00000000;
    static std::uint64_t const rpc_bx_window_mask_ = 0x000003e000000000;
    static std::uint64_t const ho_bx_window_mask_  = 0x00001c0000000000;

    static unsigned int const dt_bx_window_offset_  = 32;
    static unsigned int const rpc_bx_window_offset_ = 37;
    static unsigned int const ho_bx_window_offset_  = 42;

public:
    Header();
    Header(std::uint64_t const record[2]);
    Header(unsigned int amc_number
           , unsigned int event_counter
           , unsigned int bx_counter
           , unsigned int data_length
           , unsigned int orbit_counter
           , unsigned int board_id
           , unsigned int dt_bx_window = 0x1f
           , int rpc_bx_min = 10, int rpc_bx_max = 5 // values for !hasRPCBXWindows
           , unsigned int ho_bx_window = 0x7);
    ~Header();

    void reset();

    bool hasDTBXWindow() const;
    unsigned int getDTBXWindow() const;
    bool hasRPCBXWindow() const;
    int getRPCBXMin() const;
    int getRPCBXMax() const;
    bool hasHOBXWindow() const;
    unsigned int getHOBXWindow() const;

    void setDTBXWindow(unsigned int bx_window = 0x1f);
    void setRPCBXWindow(int bx_min = 10, int bx_max = 5); // values for !hasRPCBXWindows
    void setHOBXWindow(unsigned int bx_window = 0x7);
};

class Record
{
public:
    /** Some of the types as defined in TwinMux_uROS_payload_v12.xlsx  **/
    static unsigned int const rpc_first_type_  = 0;
    static unsigned int const rpc_second_type_ = 1;
    static unsigned int const error_type_      = 2;
    static unsigned int const unknown_type_    = 3;

    static std::uint64_t const rpc_first_identifier_mask_      = 0xf000000000000000;
    static std::uint64_t const rpc_first_identifier_           = 0x9000000000000000;
    static std::uint64_t const rpc_second_identifier_mask_     = 0xf000000000000000;
    static std::uint64_t const rpc_second_identifier_          = 0xe000000000000000;
    static std::uint64_t const error_identifier_mask_          = 0xf000000000000000;
    static std::uint64_t const error_identifier_               = 0xf000000000000000;

public:
    Record(std::uint64_t const record = 0x0);

    static unsigned int getType(std::uint64_t const record);
    unsigned int getType() const;

    void set(std::uint64_t const record);
    void reset();

    std::uint64_t const & getRecord() const;

protected:
    std::uint64_t record_;
};

class RPCLinkRecord
{
protected:
    static std::uint32_t const da_mask_             = 0x080000;
    static std::uint32_t const de_mask_             = 0x040000;
    static std::uint32_t const eod_mask_            = 0x020000;
    static std::uint32_t const delay_mask_          = 0x01c000;
    static std::uint32_t const link_board_mask_     = 0x003000;
    static std::uint32_t const connector_mask_      = 0x000e00;
    static std::uint32_t const partition_mask_      = 0x000100;
    static std::uint32_t const partition_data_mask_ = 0x0000ff;

    static unsigned int const delay_offset_          = 14;
    static unsigned int const link_board_offset_     = 12;
    static unsigned int const connector_offset_      = 9;
    static unsigned int const partition_offset_      = 8;
    static unsigned int const partition_data_offset_ = 0;

    friend class RPCRecord;

public:
    RPCLinkRecord(std::uint32_t const record = da_mask_);

    void set(std::uint32_t const record = da_mask_);
    void reset();

    std::uint32_t const & getRecord() const;

    bool isAcknowledge() const;
    bool isError() const;
    bool isEOD() const;
    unsigned int getDelay() const;
    unsigned int getLinkBoard() const;
    unsigned int getConnector() const;
    unsigned int getPartition() const;
    std::uint8_t getPartitionData() const;

    void setAcknowledge(bool da);
    void setError(bool de);
    void setEOD(bool eod);
    void setDelay(unsigned int delay);
    void setLinkBoard(unsigned int link_board);
    void setConnector(unsigned int connector);
    void setPartition(unsigned int partition);
    void setPartitionData(std::uint8_t data);

protected:
    std::uint32_t record_;
};

class RPCBXRecord
{
protected:
    static std::uint8_t const bc0_mask_ = 0x04;
    static std::uint8_t const bcn_mask_ = 0x03;
    static unsigned int const bcn_offset_ = 0;

public:
    RPCBXRecord(std::uint8_t const record = 0x0);

    void set(std::uint8_t const record = 0x0);
    void reset();

    std::uint8_t const & getRecord() const;

    bool isBC0() const;
    unsigned int getBXCounter() const;

    void setBC0(bool bc0);
    void setBXCounter(unsigned int bcn);

protected:
    std::uint8_t record_;
};

class RPCRecord
{
protected:
    static std::uint64_t const bx_offset_mask_ = 0x0f00000000000000;
    static std::uint64_t const overflow_mask_  = 0x0080000000000000;

    static unsigned int const bx_offset_offset_ = 56;

    static std::uint64_t const link_record_mask_ = 0x0fffff;
    static unsigned int const link_record_word_[5];
    static unsigned int const link_record_offset_[5];

    static std::uint64_t const bx_record_mask_ = 0x07;
    static unsigned int const bx_record_offset_[5];

public:
    RPCRecord();
    RPCRecord(std::uint64_t const record[2]);

    void set(unsigned int word, std::uint64_t const record);
    void reset();

    std::uint64_t const * getRecord() const;

    int getBXOffset() const;
    bool hasOverflow() const;
    RPCBXRecord getRPCBXRecord(unsigned int link) const;
    RPCLinkRecord getRPCLinkRecord(unsigned int link) const;

    void setBXOffset(int bx_offset);
    void setOverflow (bool overflow);
    void setRPCBXRecord(unsigned int link, RPCBXRecord const & bx_record);
    void setRPCLinkRecord(unsigned int link, RPCLinkRecord const & link_record);

protected:
    std::uint64_t record_[2];
};

} // namespace rpctwinmux

#include "EventFilter/RPCRawToDigi/interface/RPCTwinMuxRecord.icc"

#endif // EventFilter_RPCRawToDigi_RPCTwinMuxRecord_h
