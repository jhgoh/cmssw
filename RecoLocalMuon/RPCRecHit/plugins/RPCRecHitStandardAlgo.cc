/*
 *  See header file for a description of this class.
 *
 *  \author M. Maggi -- INFN
 */

#include "RPCCluster.h"
#include "RPCRecHitStandardAlgo.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/Exception.h"

// First Step
bool RPCRecHitStandardAlgo::compute(const RPCRoll& roll,
                                    const RPCCluster& cluster,
                                    LocalPoint& Point,
                                    LocalError& error,
                                    float& time,
                                    float& timeErr) const {
  // Cluster centre in the measurement frame: (strip, y/stripLength)
  const auto& topo = roll.specificTopology();
  const float stripLen = topo.stripLength();
  const float strip = (cluster.firstStrip() + cluster.lastStrip()) / 2.f - 0.5f;
  const float cls = cluster.clusterSize();

  // Across-strip error: flat over the cluster width.
  // Along-strip error: flat over the strip length, or with y the maximum estimate
  // (distance to the boundary)/sqrt(3) which is consistent with the former at y=0.
  float y = 0, vv = 1.f / 12;
  if (cluster.hasY()) {
    y = cluster.y();
    const float maxDy = stripLen / 2 - std::abs(y);
    vv = maxDy * maxDy / (3 * stripLen * stripLen);
  }

  // The topology takes care of the strip tilt, the local pitch and the x-y correlation
  const MeasurementPoint mp(strip, y / stripLen);
  Point = topo.localPosition(mp);
  error = topo.localError(mp, MeasurementError(cls * cls / 12, 0, vv));

  if (cluster.hasTime()) {
    time = cluster.time();
    timeErr = cluster.timeRMS();
  } else {
    time = 0;
    timeErr = -1;
  }

  return true;
}

bool RPCRecHitStandardAlgo::compute(const RPCRoll& roll,
                                    const RPCCluster& cl,
                                    const float& angle,
                                    const GlobalPoint& globPos,
                                    LocalPoint& Point,
                                    LocalError& error,
                                    float& time,
                                    float& timeErr) const {
  this->compute(roll, cl, Point, error, time, timeErr);
  return true;
}
