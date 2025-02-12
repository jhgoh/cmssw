#ifndef Validation_RPCRecHits_RPCValidHistograms_H
#define Validation_RPCRecHits_RPCValidHistograms_H

#include "DQMServices/Core/interface/DQMStore.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <string>

struct RPCValidHistograms {
  typedef dqm::legacy::DQMStore DQMStore;
  typedef dqm::legacy::MonitorElement MonitorElement;

  typedef MonitorElement *MEP;

  RPCValidHistograms() { booked_ = false; };

  void bookHistograms(DQMStore::IBooker &booker, const std::string &subDir);

  // Hit properties
  MEP clusterSize, clusterSizeBarrel, clusterSizeEndcap;
  MEP avgClusterSize, avgClusterSizeBarrel, avgClusterSizeEndcap;

  MEP nRefHitBarrel, nRefHitEndcap;
  MEP nRecHitBarrel, nRecHitEndcap;
  MEP nMatchHitBarrel, nMatchHitEndcap;

  MEP timeBarrel, timeEndcap, timeIRPC, timeCRPC;

  // Occupancy 1D
  MEP refHitOccupancyBarrel_wheel, refHitOccupancyEndcap_disk, refHitOccupancyBarrel_station;
  MEP recHitOccupancyBarrel_wheel, recHitOccupancyEndcap_disk, recHitOccupancyBarrel_station;
  MEP matchOccupancyBarrel_wheel, matchOccupancyEndcap_disk, matchOccupancyBarrel_station;

  // Occupancy 2D
  MEP refHitOccupancyBarrel_wheel_station, refHitOccupancyEndcap_disk_ring;
  MEP recHitOccupancyBarrel_wheel_station, recHitOccupancyEndcap_disk_ring;
  MEP matchOccupancyBarrel_wheel_station, matchOccupancyEndcap_disk_ring;

  // Residuals
  MEP resBarrel, resEndcap;
  MEP resBarrel_W0, resBarrel_W1, resBarrel_W2;
  MEP resBarrel_S1L1, resBarrel_S1L2, resBarrel_S2L1, resBarrel_S2L2, resBarrel_S3, resBarrel_S4;
  MEP resEndcap_D1, resEndcap_D2, resEndcap_D3, resEndcap_D4;
  MEP resEndcap_R1, resEndcap_R2, resEndcap_R3;

  // Pulls
  MEP pullBarrel, pullEndcap;
  MEP pullBarrel_W0, pullBarrel_W1, pullBarrel_W2;
  MEP pullBarrel_S1L1, pullBarrel_S1L2, pullBarrel_S2L1, pullBarrel_S2L2, pullBarrel_S3, pullBarrel_S4;
  MEP pullEndcap_D1, pullEndcap_D2, pullEndcap_D3, pullEndcap_D4;
  MEP pullEndcap_R1, pullEndcap_R2, pullEndcap_R3;

private:
  bool booked_;
};

#endif
