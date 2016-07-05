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

DTSegtoRPC::DTSegtoRPC(const DTRecSegment4DCollection * all4DSegments, const edm::EventSetup& iSetup, bool debug,double eyr){

  //By now hard coded parameters
  const double MinCosAng=0.85;
  const double MaxD=80.;
  const double MaxDrb4=150.;
  const double MaxDistanceBetweenSegments = 150;

  _ThePoints = std::make_unique<RPCRecHitCollection>();

  if ( all4DSegments->size() > 8 ) {
    if(debug) std::cout<<"Too many segments in this event we are not doing the extrapolation"<<std::endl;
    return;
  }

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

    if(debug) std::cout<<"DT  \t \t This Segment is in Chamber id: "<<DTId<<std::endl;
    if(debug) std::cout<<"DT  \t \t Number of segments in this DT = "<<DTSegmentCounter[DTId]<<std::endl;
    if(debug) std::cout<<"DT  \t \t Is the only one in this DT? and is not in the 4th Station?"<<std::endl;

    if ( DTSegmentCounter[DTId]!=1 or DTId.station()==4 ) {
      if(debug) std::cout<<"DT \t \t More than one segment in this chamber, or we are in Station 4"<<std::endl;
      continue;
    }

    int dtWheel = DTId.wheel();
    int dtStation = DTId.station();
    int dtSector = DTId.sector();

    LocalPoint segmentPosition= segment.localPosition();
    LocalVector segmentDirection=segment.localDirection();

    const GeomDet* gdet=dtGeo->idToDet(segment.geographicalId());
    const BoundPlane & DTSurface = gdet->surface();

    //check if the dimension of the segment is 4

    if(debug) std::cout<<"DT  \t \t Is the segment 4D?"<<std::endl;

    if(segment.dimension()!=4){
      if(debug) std::cout<<"DT  \t \t no"<<std::endl;
      continue;
    }

    if(debug) std::cout<<"DT  \t \t yes"<<std::endl;
    if(debug) std::cout<<"DT  \t \t DT Segment Dimension "<<segment.dimension()<<std::endl;

    float Xo=segmentPosition.x();
    float Yo=segmentPosition.y();
    float Zo=segmentPosition.z();
    float dx=segmentDirection.x();
    float dy=segmentDirection.y();
    float dz=segmentDirection.z();

    if(debug) std::cout<<"Creating the DTIndex"<<std::endl;
    DTStationIndex theindex(0,dtWheel,dtSector,dtStation);
    if(debug) std::cout<<"Getting the Rolls for the given index"<<std::endl;
    std::set<RPCDetId> rollsForThisDT = dtMap->getRolls(theindex);

    if(debug) std::cout<<"DT  \t \t Number of rolls for this DT = "<<rollsForThisDT.size()<<std::endl;

    assert(rollsForThisDT.size()>=1);

    if(debug) std::cout<<"DT  \t \t Loop over all the rolls asociated to this DT"<<std::endl;
    for ( auto iteraRoll : rollsForThisDT ) {
      const RPCRoll* rollasociated = rpcGeo->roll(iteraRoll);
      const RPCDetId rpcId = rollasociated->id();
      const BoundPlane & RPCSurface = rollasociated->surface();

      RPCGeomServ rpcsrv(rpcId);
      const std::string nameRoll = rpcsrv.name();

      if(debug) std::cout<<"DT  \t \t \t RollName: "<<nameRoll<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t Doing the extrapolation to this roll"<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t DT Segment Direction in DTLocal "<<segmentDirection<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t DT Segment Point in DTLocal "<<segmentPosition<<std::endl;

      const GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));

      const LocalPoint CenterRollinDTFrame = DTSurface.toLocal(CenterPointRollGlobal);

      if(debug) std::cout<<"DT  \t \t \t Center (0,0,0) Roll In DTLocal"<<CenterRollinDTFrame<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t Center (0,0,0) of the Roll in Global"<<CenterPointRollGlobal<<std::endl;

      float D=CenterRollinDTFrame.z();

      float X=Xo+dx*D/dz;
      float Y=Yo+dy*D/dz;
      float Z=D;

      const RectangularStripTopology* topo= dynamic_cast<const RectangularStripTopology*> (&(rollasociated->topology()));
      LocalPoint xmin = topo->localPosition(0.);
      if(debug) std::cout<<"DT  \t \t \t xmin of this  Roll "<<xmin<<"cm"<<std::endl;
      LocalPoint xmax = topo->localPosition((float)rollasociated->nstrips());
      if(debug) std::cout<<"DT  \t \t \t xmax of this  Roll "<<xmax<<"cm"<<std::endl;
      float rsize = fabs( xmax.x()-xmin.x() );
      if(debug) std::cout<<"DT  \t \t \t Roll Size "<<rsize<<"cm"<<std::endl;
      float stripl = topo->stripLength();

      float stripw = topo->pitch();

      if(debug) std::cout<<"DT  \t \t \t Strip Lenght "<<stripl<<"cm"<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t Strip Width "<<stripw<<"cm"<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t X Predicted in DTLocal= "<<X<<"cm"<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t Y Predicted in DTLocal= "<<Y<<"cm"<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t Z Predicted in DTLocal= "<<Z<<"cm"<<std::endl;

      float extrapolatedDistance = sqrt((X-Xo)*(X-Xo)+(Y-Yo)*(Y-Yo)+(Z-Zo)*(Z-Zo));

      if(debug) std::cout<<"DT  \t \t \t Is the distance of extrapolation less than MaxD? ="<<extrapolatedDistance<<"cm"<<"MaxD="<<MaxD<<"cm"<<std::endl;

      if ( extrapolatedDistance > MaxD ) {
        if(debug) std::cout<<"DT \t \t \t No, Exrtrapolation too long!, canceled"<<std::endl;
        continue;
      }

      if(debug) std::cout<<"DT  \t \t \t yes"<<std::endl;
      GlobalPoint GlobalPointExtrapolated = DTSurface.toGlobal(LocalPoint(X,Y,Z));
      if(debug) std::cout<<"DT  \t \t \t Point ExtraPolated in Global"<<GlobalPointExtrapolated<< std::endl;
      LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);

      if(debug) std::cout<<"DT  \t \t \t Point Extrapolated in RPCLocal"<<PointExtrapolatedRPCFrame<< std::endl;
      if(debug) std::cout<<"DT  \t \t \t Corner of the Roll = ("<<rsize*eyr<<","<<stripl*eyr<<")"<<std::endl;
      if(debug) std::cout<<"DT \t \t \t Info About the Point Extrapolated in X Abs ("<<fabs(PointExtrapolatedRPCFrame.x())<<","
        <<fabs(PointExtrapolatedRPCFrame.y())<<","<<fabs(PointExtrapolatedRPCFrame.z())<<")"<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t Does the extrapolation go inside this roll?"<<std::endl;

      if ( fabs(PointExtrapolatedRPCFrame.z()) >= 1. or
           fabs(PointExtrapolatedRPCFrame.x()) >= rsize*eyr or
           fabs(PointExtrapolatedRPCFrame.y()) >= stripl*eyr ) {
        if(debug) std::cout<<"DT \t \t \t \t No the prediction is outside of this roll"<<std::endl;
        continue;
      }

      if(debug) std::cout<<"DT  \t \t \t \t yes"<<std::endl;
      if(debug) std::cout<<"DT  \t \t \t \t Creating the RecHit"<<std::endl;
      RPCRecHit RPCPoint(rpcId,0,PointExtrapolatedRPCFrame);
      if(debug) std::cout<<"DT  \t \t \t \t Clearing the vector"<<std::endl;
      RPCPointVector.clear();
      if(debug) std::cout<<"DT  \t \t \t \t Pushing back"<<std::endl;
      RPCPointVector.push_back(RPCPoint);
      if(debug) std::cout<<"DT  \t \t \t \t Putting the vector"<<std::endl;
      _ThePoints->put(rpcId,RPCPointVector.begin(),RPCPointVector.end());

      if(debug) std::cout<<"DT \t \t \t \t Filling container with "<<nameRoll
                         <<" Point.x="<<PointExtrapolatedRPCFrame.x()<<" Point.y="<<PointExtrapolatedRPCFrame.y()<<" size="<<RPCPointVector.size()<<std::endl;
    }//loop over all the rolls asociated
  }

  if ( all4DSegments->size() == 0 ) {
    if(debug) std::cout<<"MB4 \t This event doesn't have 4D Segment"<<std::endl;
    return;
  }
  if(debug) std::cout<<"MB4 \t \t Loop Over all4DSegments "<<all4DSegments->size()<<std::endl;
  extrapolatedRolls.clear();
  for (auto& segment : *all4DSegments ) {

    DTChamberId DTId = segment.chamberId();

    if(debug) std::cout<<"MB4 \t \t This Segment is in Chamber id: "<<DTId<<std::endl;
    if(debug) std::cout<<"MB4 \t \t Number of segments in this DT = "<<DTSegmentCounter[DTId]<<std::endl;
    if(debug) std::cout<<"MB4 \t \t \t Is the only one in this DT? and is in the Station 4?"<<std::endl;

    if ( DTSegmentCounter[DTId] != 1 or DTId.station()!=4 ){
      if(debug) std::cout<<"MB4 \t \t \t \t There is not just one segment or is not in station 4"<<std::endl;
      continue;
    }

    if(debug) std::cout<<"MB4 \t \t \t yes"<<std::endl;
    int dtWheel = DTId.wheel();
    int dtStation = DTId.station();
    int dtSector = DTId.sector();

    LocalPoint segmentPosition= segment.localPosition();
    LocalVector segmentDirection=segment.localDirection();

    if(debug) std::cout<<"MB4 \t \t \t \t The Segment in MB4 is 2D?"<<std::endl;
    if(segment.dimension()!=2){
      if(debug) std::cout<<"MB4 \t \t \t Is NOT a 2D Segment"<<std::endl;
      continue;
    }
    if(debug) std::cout<<"MB4 \t \t \t \t yes"<<std::endl;
    LocalVector segmentDirectionMB4=segmentDirection;
    LocalPoint segmentPositionMB4=segmentPosition;

    const BoundPlane& DTSurface4 = dtGeo->idToDet(DTId)->surface();

    if(debug) std::cout<<"MB4 \t \t \t \t Loop on segments in =sector && MB3 && adjacent sectors && y dim=4"<<std::endl;
    for ( auto& segMB3 : *all4DSegments ) {

      DTChamberId dtid3 = segMB3.chamberId();

      if(debug) std::cout<<"MB4  \t \t \t \t Segment in Chamber ="<<dtid3<<std::endl;

      if ( distsector(dtid3.sector(),DTId.sector()) > 1 //The DT sector could be 13 or 14 and because is corrected in the calculation of the distance.
           || distwheel(dtid3.wheel(),DTId.wheel()) > 1 //The we could have segments in neighbohr wheels in pp collisions
           || dtid3.station()!=3
           || DTSegmentCounter[dtid3] != 1
           || segMB3.dimension()!=4){
        if(debug) std::cout<<"MB4 \t \t \t No the same station or same wheel or segment dim in mb3 not 4D"<<std::endl;
        continue;
      }

      if(debug) std::cout<<"MB4  \t \t \t \t distsector ="<<distsector(dtid3.sector(),DTId.sector())<<std::endl;
      if(debug) std::cout<<"MB4  \t \t \t \t distwheel ="<<distwheel(dtid3.wheel(),DTId.wheel())<<std::endl;

      const GeomDet* gdet3=dtGeo->idToDet(segMB3.geographicalId());
      const BoundPlane & DTSurface3 = gdet3->surface();

      LocalVector segmentDirectionMB3 =  segMB3.localDirection();
      GlobalPoint segmentPositionMB3inGlobal = DTSurface3.toGlobal(segMB3.localPosition());
      GlobalPoint segmentPositionMB4inGlobal = DTSurface4.toGlobal(segmentPosition);

      //LocalVector segDirMB4inMB3Frame=DTSurface3.toLocal(DTSurface4.toGlobal(segmentDirectionMB4));
      LocalVector segDirMB3inMB4Frame=DTSurface4.toLocal(DTSurface3.toGlobal(segmentDirectionMB3));

      GlobalVector segDirMB4inGlobalFrame=DTSurface4.toGlobal(segmentDirectionMB4);
      GlobalVector segDirMB3inGlobalFrame=DTSurface3.toGlobal(segmentDirectionMB3);

      float dx=segDirMB4inGlobalFrame.x();
      float dy=segDirMB4inGlobalFrame.y();

      float dx3=segDirMB3inGlobalFrame.x();
      float dy3=segDirMB3inGlobalFrame.y();

      double cosAng=fabs(dx*dx3+dy*dy3/sqrt((dx3*dx3+dy3*dy3)*(dx*dx+dy*dy)));

      if(debug) std::cout<<"MB4 \t \t \t \t cosAng"<<cosAng<<"Beetween "<<dtid3<<" and "<<DTId<<std::endl;

      if(debug){
        std::cout<<"MB4 \t \t \t \t dx="<<dx<<" dy="<<dy<<std::endl;
        std::cout<<"MB4 \t \t \t \t dx3="<<dx3<<" dy3="<<dy<<std::endl;
        std::cout<<"MB4 \t \t \t \t cosAng="<<cosAng<<std::endl;
      }

      float DistanceBetweenSegments = ((segmentPositionMB3inGlobal) - (segmentPositionMB4inGlobal)).mag();

      if ( cosAng<=MinCosAng or DistanceBetweenSegments >= MaxDistanceBetweenSegments){
        if(debug) std::cout<<"MB4 \t \t \t \t I found segments in MB4 and MB3 adjacent wheel and/or sector but not compatibles, Diferent Directions"<<std::endl;
        continue;
      }

      if(debug) std::cout<<"MB4 \t \t \t \t Distance between segments="<<DistanceBetweenSegments<<std::endl;

      if(debug) std::cout<<"MB4 \t \t We found compatible Segments (similar direction and close enough) in "<<dtid3<<" and "<<DTId<<std::endl;

      if      ( dtSector==13 ) dtSector=4;
      else if ( dtSector==14 ) dtSector=10;

      if(debug) std::cout<<"Creating the DTIndex"<<std::endl;
      DTStationIndex theindex(0,dtWheel,dtSector,dtStation);
      if(debug) std::cout<<"Getting the Rolls for the given index"<<std::endl;
      std::set<RPCDetId> rollsForThisDT = dtMap->getRolls(theindex);

      if(debug) std::cout<<"MB4 \t \t Number of rolls for this DT = "<<rollsForThisDT.size()<<std::endl;

      assert(rollsForThisDT.size()>=1);

      if(debug) std::cout<<"MB4  \t \t Loop over all the rolls asociated to this DT"<<std::endl;
      for (auto iteraRoll : rollsForThisDT ) {
        const RPCRoll* rollasociated = rpcGeo->roll(iteraRoll); //roll asociado a MB4
        RPCDetId rpcId = rollasociated->id();
        const BoundPlane & RPCSurfaceRB4 = rollasociated->surface(); //surface MB4

        RPCGeomServ rpcsrv(rpcId);
        std::string nameRoll = rpcsrv.name();

        if(debug) std::cout<<"MB4  \t \t \t RollName: "<<nameRoll<<std::endl;
        if(debug) std::cout<<"MB4  \t \t \t Doing the extrapolation to this roll"<<std::endl;

        GlobalPoint CenterPointRollGlobal=RPCSurfaceRB4.toGlobal(LocalPoint(0,0,0));
        LocalPoint CenterRollinMB4Frame = DTSurface4.toLocal(CenterPointRollGlobal); //In MB4
        LocalPoint segmentPositionMB3inMB4Frame = DTSurface4.toLocal(segmentPositionMB3inGlobal); //In MB4
        //LocalPoint segmentPositionMB3inRB4Frame = RPCSurfaceRB4.toLocal(segmentPositionMB3inGlobal); //In MB4
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

        if(debug) std::cout<<"MB4 \t \t \t The distance to extrapolate in Y from MB3 is "<<Dy<<"cm"<<std::endl;

        float Y=Yo34+dy34*Dy/dz34;//In MB4 Frame

        const auto topo = dynamic_cast<const RectangularStripTopology*>(&(rollasociated->topology())); //Topology roll asociated MB4
        LocalPoint xmin = topo->localPosition(0.);
        LocalPoint xmax = topo->localPosition((float)rollasociated->nstrips());
        float rsize = fabs( xmax.x()-xmin.x() );
        float stripl = topo->stripLength();
        float stripw = topo->pitch();


        if(debug) std::cout<<"MB4 \t \t \t Strip Lenght "<<stripl<<"cm"<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t Strip Width "<<stripw<<"cm"<<std::endl;

        if(debug) std::cout<<"MB4 \t \t \t X Predicted in MB4DTLocal= "<<X<<"cm"<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t Y Predicted in MB4DTLocal= "<<Y<<"cm"<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t Z Predicted in MB4DTLocal= "<<Z<<"cm"<<std::endl;

        float extrapolatedDistance = sqrt((Y-Yo34)*(Y-Yo34)+Dy*Dy);

        if(debug) std::cout<<"MB4 \t \t \t segmentPositionMB3inMB4Frame"<<segmentPositionMB3inMB4Frame<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t segmentPositionMB4inMB4Frame"<<segmentPosition<<std::endl;

        if(debug) std::cout<<"MB4 \t \t \t segmentDirMB3inMB4Frame"<<segDirMB3inMB4Frame<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t segmentDirMB4inMB4Frame"<<segmentDirectionMB4<<std::endl;

        if(debug) std::cout<<"MB4 \t \t \t CenterRB4PositioninMB4Frame"<<CenterRollinMB4Frame<<std::endl;

        if(debug) std::cout<<"MB4 \t \t \t Is the extrapolation distance ="<<extrapolatedDistance<<"less than "<<MaxDrb4<<std::endl;


        if ( extrapolatedDistance>MaxDrb4 ) {
          if(debug) std::cout<<"MB4 \t \t \t No, Exrtrapolation too long!, canceled"<<std::endl;
          continue;
        }
        if(debug) std::cout<<"MB4 \t \t \t yes"<<std::endl;

        GlobalPoint GlobalPointExtrapolated = DTSurface4.toGlobal(LocalPoint(X,Y,Z));

        if(debug) std::cout<<"MB4 \t \t \t Point ExtraPolated in Global"<<GlobalPointExtrapolated<< std::endl;

        LocalPoint PointExtrapolatedRPCFrame = RPCSurfaceRB4.toLocal(GlobalPointExtrapolated);

        if(debug) std::cout<<"MB4 \t \t \t Point Extrapolated in RPCLocal"<<PointExtrapolatedRPCFrame<< std::endl;
        if(debug) std::cout<<"MB4 \t \t \t Corner of the Roll = ("<<rsize*eyr<<","<<stripl*eyr<<")"<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t Info About the Point Extrapolated in X Abs ("<<fabs(PointExtrapolatedRPCFrame.x())<<","
                           <<fabs(PointExtrapolatedRPCFrame.y())<<","<<fabs(PointExtrapolatedRPCFrame.z())<<")"<<std::endl;

        if(debug) std::cout<<"MB4 \t \t \t Does the extrapolation go inside this roll?"<<std::endl;

        if ( fabs(PointExtrapolatedRPCFrame.z()) >= 5. or
            fabs(PointExtrapolatedRPCFrame.x()) >= rsize*eyr or
            fabs(PointExtrapolatedRPCFrame.y()) >= stripl*eyr ) {
          if(debug) std::cout<<"MB4 \t \t \t \t No the prediction is outside of this roll"<<std::endl;
          continue;
        }
        if(debug) std::cout<<"MB4 \t \t \t \t yes"<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t \t Creating the RecHit"<<std::endl;
        RPCRecHit RPCPointMB4(rpcId,0,PointExtrapolatedRPCFrame);
        if(debug) std::cout<<"MB4 \t \t \t \t Clearing the RPCPointVector"<<std::endl;
        RPCPointVector.clear();
        if(debug) std::cout<<"MB4 \t \t \t \t Pushing Back"<<std::endl;
        RPCPointVector.push_back(RPCPointMB4);
        if(debug) std::cout<<"MB4 \t \t \t \t Putting for "<<rpcId<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t \t Filling container with "<<nameRoll
                           <<" Point.x="<<PointExtrapolatedRPCFrame.x()<<" Point.y="<<PointExtrapolatedRPCFrame.y()<<" size="<<RPCPointVector.size()<<std::endl;
        if(debug) std::cout<<"MB4 \t \t \t \t Number of rolls already extrapolated in RB4 = "<<extrapolatedRolls.size()<<std::endl;
        if(find (extrapolatedRolls.begin(),extrapolatedRolls.end(),rpcId.rawId()) != extrapolatedRolls.end()){
          if(debug) std::cout<<"MB4 \t \t \t \t roll already extrapolated "<<rpcId<<std::endl;
          continue;
        }

        extrapolatedRolls.push_back(rpcId.rawId());
        _ThePoints->put(rpcId,RPCPointVector.begin(),RPCPointVector.end());
        if(debug) std::cout<<"MB4 \t \t \t \t Extrapolations done after this point = "<<extrapolatedRolls.size()<<std::endl;
        if(debug) for(uint32_t m=0;m<extrapolatedRolls.size();m++) std::cout<<"MB4 \t \t \t \t"<< extrapolatedRolls.at(m)<<std::endl;
      }
    }
  }
}

DTSegtoRPC::~DTSegtoRPC() { }
