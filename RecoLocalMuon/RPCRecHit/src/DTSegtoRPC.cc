#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "RecoLocalMuon/RPCRecHit/interface/DTSegtoRPC.h"
#include "RecoLocalMuon/RPCRecHit/src/DTObjectMap.h"
#include "RecoLocalMuon/RPCRecHit/src/DTStationIndex.h"

#include <ctime>

int distsector(int sector1,int sector2){
  if(sector1==13) sector1=4;
  if(sector1==14) sector1=10;

  if(sector2==13) sector2=4;
  if(sector2==14) sector2=10;

  int distance = std::abs(sector1 - sector2);
  if(distance>6) distance = 12-distance;
  return distance;
}

int distwheel(int wheel1,int wheel2){
  int distance = std::abs(wheel1 - wheel2);
  return distance;
}

DTSegtoRPC::DTSegtoRPC(const DTRecSegment4DCollection * all4DSegments, const edm::EventSetup& iSetup, double eyr){

  //By now hard coded parameters
  const double MinCosAng=0.85;
  const double MaxD=80.;
  const double MaxDrb4=160.;
  const double MaxDistanceBetweenSegments = 160;

  _ThePoints.reset(new RPCRecHitCollection());

  if ( all4DSegments->size() > 8 ) return;

  edm::ESHandle<RPCGeometry> rpcGeo;
  edm::ESHandle<DTGeometry> dtGeo;
  edm::ESHandle<DTObjectMap> dtMap;

  iSetup.get<MuonGeometryRecord>().get(rpcGeo);
  iSetup.get<MuonGeometryRecord>().get(dtGeo);
  iSetup.get<MuonGeometryRecord>().get(dtMap);

  std::map<DTChamberId,int> DTSegmentCounter;
  for ( auto& segment : *all4DSegments ) {
    ++DTSegmentCounter[segment.chamberId()];
  }

  for ( auto& segment : *all4DSegments ) {

    const DTChamberId DTId = segment.chamberId();

    if ( DTSegmentCounter[DTId]!=1 or DTId.station()==4 ) continue;

    int dtWheel = DTId.wheel();
    int dtStation = DTId.station();
    int dtSector = DTId.sector();

    LocalPoint segmentPosition= segment.localPosition();
    LocalVector segmentDirection=segment.localDirection();

    const GeomDet* gdet=dtGeo->idToDet(segment.geographicalId());
    const BoundPlane& DTSurface = gdet->surface();

    //check if the dimension of the segment is 4

    if(segment.dimension()!=4) continue;

    float Xo=segmentPosition.x();
    float Yo=segmentPosition.y();
    float Zo=segmentPosition.z();
    float dx=segmentDirection.x();
    float dy=segmentDirection.y();
    float dz=segmentDirection.z();

    DTStationIndex theindex(0,dtWheel,dtSector,dtStation);
    std::set<RPCDetId> rollsForThisDT = dtMap->getRolls(theindex);

    assert(rollsForThisDT.size()>=1);

    for ( auto iteraRoll : rollsForThisDT ) {
      const RPCRoll* rollasociated = rpcGeo->roll(iteraRoll);
      const RPCDetId rpcId = rollasociated->id();
      const BoundPlane & RPCSurface = rollasociated->surface();

      const GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));
      const LocalPoint CenterRollinDTFrame = DTSurface.toLocal(CenterPointRollGlobal);

      float D=CenterRollinDTFrame.z();

      float X=Xo+dx*D/dz;
      float Y=Yo+dy*D/dz;
      float Z=D;

      auto topo= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
      LocalPoint xmin = topo->localPosition(0.);
      LocalPoint xmax = topo->localPosition((float)rollasociated->nstrips());
      float rsize = std::abs( xmax.x()-xmin.x() );
      float stripl = topo->stripLength();

      float extrapolatedDistance = sqrt((X-Xo)*(X-Xo)+(Y-Yo)*(Y-Yo)+(Z-Zo)*(Z-Zo));

      if ( extrapolatedDistance > MaxD ) continue;

      GlobalPoint GlobalPointExtrapolated = DTSurface.toGlobal(LocalPoint(X,Y,Z));
      LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);

      if ( std::abs(PointExtrapolatedRPCFrame.z()) >= 1. or
           std::abs(PointExtrapolatedRPCFrame.x()) >= rsize*eyr or
           std::abs(PointExtrapolatedRPCFrame.y()) >= stripl*eyr ) {
        continue;
      }

      RPCRecHit RPCPoint(rpcId,0,PointExtrapolatedRPCFrame);
      RPCPointVector.clear();
      RPCPointVector.push_back(RPCPoint);
      _ThePoints->put(rpcId,RPCPointVector.begin(),RPCPointVector.end());

    }//loop over all the rolls asociated
  }

  if ( all4DSegments->size() == 0 ) return;
  extrapolatedRolls.clear();
  for (auto& segment : *all4DSegments ) {

    DTChamberId DTId = segment.chamberId();

    if ( DTSegmentCounter[DTId] != 1 or DTId.station()!=4 ) continue;

    int dtWheel = DTId.wheel();
    int dtStation = DTId.station();
    int dtSector = DTId.sector();

    LocalPoint segmentPosition= segment.localPosition();
    LocalVector segmentDirection=segment.localDirection();

    if(segment.dimension()!=2) continue;
    LocalVector segmentDirectionMB4=segmentDirection;
    LocalPoint segmentPositionMB4=segmentPosition;

    const BoundPlane& DTSurface4 = dtGeo->idToDet(DTId)->surface();

    for ( auto& segMB3 : *all4DSegments ) {

      DTChamberId dtid3 = segMB3.chamberId();

      if ( distsector(dtid3.sector(),DTId.sector()) > 1 //The DT sector could be 13 or 14 and because is corrected in the calculation of the distance.
           || distwheel(dtid3.wheel(),DTId.wheel()) > 1 //The we could have segments in neighbohr wheels in pp collisions
           || dtid3.station()!=3
           || DTSegmentCounter[dtid3] != 1
           || segMB3.dimension()!=4){
        continue;
      }

      const GeomDet* gdet3=dtGeo->idToDet(segMB3.geographicalId());
      const BoundPlane & DTSurface3 = gdet3->surface();

      LocalVector segmentDirectionMB3 =  segMB3.localDirection();
      GlobalPoint segmentPositionMB3inGlobal = DTSurface3.toGlobal(segMB3.localPosition());
      GlobalPoint segmentPositionMB4inGlobal = DTSurface4.toGlobal(segmentPosition);

      GlobalVector segDirMB4inGlobalFrame=DTSurface4.toGlobal(segmentDirectionMB4);
      GlobalVector segDirMB3inGlobalFrame=DTSurface3.toGlobal(segmentDirectionMB3);

      float dx=segDirMB4inGlobalFrame.x();
      float dy=segDirMB4inGlobalFrame.y();

      float dx3=segDirMB3inGlobalFrame.x();
      float dy3=segDirMB3inGlobalFrame.y();

      double cosAng=std::abs(dx*dx3+dy*dy3/sqrt((dx3*dx3+dy3*dy3)*(dx*dx+dy*dy)));

      float DistanceBetweenSegments = ((segmentPositionMB3inGlobal) - (segmentPositionMB4inGlobal)).mag();
      if ( cosAng<=MinCosAng or DistanceBetweenSegments >= MaxDistanceBetweenSegments) continue;

      if      ( dtSector==13 ) dtSector=4;
      else if ( dtSector==14 ) dtSector=10;

      DTStationIndex theindex(0,dtWheel,dtSector,dtStation);
      std::set<RPCDetId> rollsForThisDT = dtMap->getRolls(theindex);

      assert(rollsForThisDT.size()>=1);

      for (auto iteraRoll : rollsForThisDT ) {
        const RPCRoll* rollasociated = rpcGeo->roll(iteraRoll); //roll asociado a MB4
        RPCDetId rpcId = rollasociated->id();
        const BoundPlane & RPCSurfaceRB4 = rollasociated->surface(); //surface MB4

        RPCGeomServ rpcsrv(rpcId);

        GlobalPoint CenterPointRollGlobal=RPCSurfaceRB4.toGlobal(LocalPoint(0,0,0));
        LocalPoint CenterRollinMB4Frame = DTSurface4.toLocal(CenterPointRollGlobal); //In MB4
        LocalPoint segmentPositionMB3inMB4Frame = DTSurface4.toLocal(segmentPositionMB3inGlobal); //In MB4
        LocalVector segmentDirectionMB3inMB4Frame = DTSurface4.toLocal(segDirMB3inGlobalFrame); //In MB4

        //The exptrapolation is done in MB4 frame. for local x and z is done from MB4,
        float Dxz=CenterRollinMB4Frame.z();
        float Xo4=segmentPositionMB4.x();
        float dxl=segmentDirectionMB4.x(); //dx local for MB4 segment in MB4 Frame
        float dzl=segmentDirectionMB4.z(); //dx local for MB4 segment in MB4 Frame

        float X=Xo4+dxl*Dxz/dzl; //In MB4 frame
        float Z=Dxz;//In MB4 frame

        //for local y is done from MB3
        float Yo34 = segmentPositionMB3inMB4Frame.y();
        float dy34 = segmentDirectionMB3inMB4Frame.y();
        float dz34 = segmentDirectionMB3inMB4Frame.z();
        float Dy=Dxz-(segmentPositionMB3inMB4Frame.z()); //Distance beetween the segment in MB3 and the RB4 surface

        float Y=Yo34+dy34*Dy/dz34;//In MB4 Frame

        const auto topo = dynamic_cast<const RectangularStripTopology*>(&(rollasociated->topology())); //Topology roll asociated MB4
        LocalPoint xmin = topo->localPosition(0.);
        LocalPoint xmax = topo->localPosition((float)rollasociated->nstrips());
        float rsize = std::abs( xmax.x()-xmin.x() );
        float stripl = topo->stripLength();

        float extrapolatedDistance = sqrt((Y-Yo34)*(Y-Yo34)+Dy*Dy);
        if ( extrapolatedDistance>MaxDrb4 ) continue;

        GlobalPoint GlobalPointExtrapolated = DTSurface4.toGlobal(LocalPoint(X,Y,Z));
        LocalPoint PointExtrapolatedRPCFrame = RPCSurfaceRB4.toLocal(GlobalPointExtrapolated);

        if ( std::abs(PointExtrapolatedRPCFrame.z()) >= 5. or
             std::abs(PointExtrapolatedRPCFrame.x()) >= rsize*eyr or
             std::abs(PointExtrapolatedRPCFrame.y()) >= stripl*eyr ) {
          continue;
        }
        RPCRecHit RPCPointMB4(rpcId,0,PointExtrapolatedRPCFrame);
        RPCPointVector.clear();
        RPCPointVector.push_back(RPCPointMB4);
        if(find (extrapolatedRolls.begin(),extrapolatedRolls.end(),rpcId.rawId()) != extrapolatedRolls.end()) continue;

        extrapolatedRolls.push_back(rpcId.rawId());
        _ThePoints->put(rpcId,RPCPointVector.begin(),RPCPointVector.end());
      }
    }
  }
}

DTSegtoRPC::~DTSegtoRPC() { }
