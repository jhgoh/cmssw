#include "DataFormats/RPCDigi/interface/RPCCPPFDigi.h"

#include <ostream>

RPCCPPFDigi::RPCCPPFDigi(RPCDetId const & rpc_id
                         , int bx
                         , unsigned int theta, unsigned int phi)
    : rpc_id_(rpc_id)
    , bx_(bx)
    , theta_(theta), phi_(phi)
{}

bool operator<(RPCCPPFDigi const & lhs, RPCCPPFDigi const & rhs)
{
    return (lhs.getRPCDetId().rawId() < rhs.getRPCDetId().rawId()
            || (!(rhs.getRPCDetId().rawId() < lhs.getRPCDetId().rawId())
                && (lhs.getBX() < rhs.getBX()
                    || (!(rhs.getBX() < lhs.getBX())
                        && (lhs.getTheta() < rhs.getTheta()
                            || (!(rhs.getTheta() < lhs.getTheta())
                                && lhs.getPhi() < rhs.getPhi()))))));
}

std::ostream & operator<<(std::ostream & os, RPCCPPFDigi const & digi)
{
    return os << "RPCCPPFDigi(" << digi.getRPCDetId().rawId()
              << ", " << digi.getBX()
              << ", " << digi.getTheta() << ", " << digi.getPhi()
              << ")";
}
