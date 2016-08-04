#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "RecoLocalMuon/RPCRecHit/interface/CSCSegtoRPC.h"
#include "RecoLocalMuon/RPCRecHit/src/CSCStationIndex.h"
#include "RecoLocalMuon/RPCRecHit/src/CSCObjectMap.h"

CSCSegtoRPC::CSCSegtoRPC(const CSCSegmentCollection * allCSCSegments, const edm::EventSetup& iSetup, double eyr){

  edm::ESHandle<RPCGeometry> rpcGeo;
  edm::ESHandle<CSCGeometry> cscGeo;
  edm::ESHandle<CSCObjectMap> cscMap;

  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  iSetup.get<MuonGeometryRecord>().get(cscGeo);
  iSetup.get<MuonGeometryRecord>().get(cscMap);

  const double MaxD=80.;

  if(debug) std::cout<<"CSC \t Number of CSC Segments in this event = "<<allCSCSegments->size()<<std::endl;
  _ThePoints = std::make_unique<RPCRecHitCollection>();

  if(allCSCSegments->size()==0) return;
  std::map<CSCDetId,int> CSCSegmentsCounter;

  int segmentsInThisEventInTheEndcap=0;
  for ( auto segment : *allCSCSegments ) {
    ++CSCSegmentsCounter[segment.cscDetId()];
    ++segmentsInThisEventInTheEndcap;
  }

  for ( auto segment : *allCSCSegments ) {
    const CSCDetId CSCId = segment.cscDetId();

    if ( CSCSegmentsCounter[CSCId]!=1 or CSCId.ring() ==1 or allCSCSegments->size()<2 ) {
      continue;
    }

    const int cscEndCap = CSCId.endcap();
    const int cscStation = CSCId.station();
    const int cscRing = CSCId.ring();
    const int rpcRegion = cscEndCap == 2 ? -1 : 1;
    const int rpcRing = cscRing==4 ? 1 : cscRing;
    const int rpcStation = cscStation;
    const int rpcSegment = CSCId.chamber();

    const LocalPoint segmentPosition= segment.localPosition();
    const LocalVector segmentDirection=segment.localDirection();
    const float dz=segmentDirection.z();

    if ( segment.dimension()!=4 or segment.nRecHits() > 10 or segment.nRecHits() < 4 ) continue;

    const float Xo=segmentPosition.x();
    const float Yo=segmentPosition.y();
    const float Zo=segmentPosition.z();
    const float dx=segmentDirection.x();
    const float dy=segmentDirection.y();

    const CSCStationIndex theindex(rpcRegion,rpcStation,rpcRing,rpcSegment);
    std::set<RPCDetId> rollsForThisCSC = cscMap->getRolls(theindex);
    const CSCChamber* TheChamber=cscGeo->chamber(CSCId);

    if(rpcRing==1) continue; //They don't exist!

    assert(rollsForThisCSC.size()>=1);

    for ( auto iteraRoll : rollsForThisCSC ) {
      const RPCRoll* rollasociated = rpcGeo->roll(iteraRoll);
      const RPCDetId rpcId = rollasociated->id();

      const BoundPlane & RPCSurface = rollasociated->surface();

      const GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));
      const LocalPoint CenterRollinCSCFrame = TheChamber->toLocal(CenterPointRollGlobal);

      const float D=CenterRollinCSCFrame.z();

      const float X=Xo+dx*D/dz;
      const float Y=Yo+dy*D/dz;
      const float Z=D;

      const auto topo = dynamic_cast<const TrapezoidalStripTopology*>(&(rollasociated->topology()));
      const LocalPoint xmin = topo->localPosition(0.);
      const LocalPoint xmax = topo->localPosition((float)rollasociated->nstrips());
      const float rsize = std::abs( xmax.x()-xmin.x() );
      const float stripl = topo->stripLength();

      const float extrapolatedDistance = sqrt((X-Xo)*(X-Xo)+(Y-Yo)*(Y-Yo)+(Z-Zo)*(Z-Zo));
      if(extrapolatedDistance < MaxD) continue;

      const GlobalPoint GlobalPointExtrapolated=TheChamber->toGlobal(LocalPoint(X,Y,Z));
      const LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);

      if ( std::abs(PointExtrapolatedRPCFrame.z()) >= 1. or
           std::abs(PointExtrapolatedRPCFrame.x()) >= rsize*eyr or
           std::abs(PointExtrapolatedRPCFrame.y()) >= stripl*eyr) continue;

      RPCRecHit RPCPoint(rpcId,0,PointExtrapolatedRPCFrame);
      RPCPointVector.clear();
      RPCPointVector.push_back(RPCPoint);
      _ThePoints->put(rpcId,RPCPointVector.begin(),RPCPointVector.end());
    }
  }
}

CSCSegtoRPC::~CSCSegtoRPC() { }
