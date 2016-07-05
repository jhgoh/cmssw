#ifndef  TRACKTORPC_H
#define  TRACKTORPC_H

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TrackRefitter/interface/TrackTransformerBase.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

#include <memory>

class TracktoRPC {
public:
  TracktoRPC(reco::TrackCollection const* alltracks, edm::EventSetup const& iSetup, bool debug, const edm::ParameterSet& iConfig, const edm::InputTag & tracklabel);
  ~TracktoRPC() = default;
  std::unique_ptr<RPCRecHitCollection> && thePoints(){ return std::move(_ThePoints); }

private:
  bool ValidRPCSurface(RPCDetId rpcid, LocalPoint LocalP, const edm::EventSetup& iSetup);

  std::unique_ptr<RPCRecHitCollection> _ThePoints;
  edm::OwnVector<RPCRecHit> RPCPointVector;
  double MaxD;

 TrackTransformerBase *theTrackTransformer;
 edm::ESHandle<Propagator> thePropagator;
};

#endif
