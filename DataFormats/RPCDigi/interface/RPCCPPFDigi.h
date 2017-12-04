#ifndef DataFormats_RPCDigi_RPCCPPFDigi_h
#define DataFormats_RPCDigi_RPCCPPFDigi_h

#include <vector>
#include <iosfwd>

#include "DataFormats/MuonDetId/interface/RPCDetId.h"

class RPCCPPFDigi
{
public:
    explicit RPCCPPFDigi(RPCDetId const & rpc_id = RPCDetId()
                         , int bx = 0
                         , unsigned int theta = 0x1f, unsigned int phi = 0x7ff);

    RPCDetId const & getRPCDetId() const;
    int getBX() const;
    unsigned int getTheta() const;
    unsigned int getPhi() const;

    bool isValid() const;

    void setRPCDetId(RPCDetId const & rpc_id);
    void setBX(int bx);
    void setTheta(int theta);
    void setPhi(int phi);

protected:
    RPCDetId rpc_id_;
    int bx_;
    unsigned int theta_, phi_;
};

bool operator<(RPCCPPFDigi const & lhs, RPCCPPFDigi const & rhs);
std::ostream & operator<<(std::ostream & os, RPCCPPFDigi const & digi);

typedef std::vector<RPCCPPFDigi> RPCCPPFDigiCollection;

#include "DataFormats/RPCDigi/interface/RPCCPPFDigi.icc"

#endif // DataFormats_RPCDigi_RPCCPPFDigi_h
