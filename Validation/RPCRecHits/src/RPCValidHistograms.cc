#include "Validation/RPCRecHits/interface/RPCValidHistograms.h"

#include "TAxis.h"

void RPCValidHistograms::bookHistograms(DQMStore::IBooker &booker, const std::string &subDir) {
  if (booked_) {
    edm::LogError("RPCValidHistograms") << "Histogram is already booked\n";
    return;
  }

  const std::string pwd = booker.pwd();
  booker.setCurrentFolder(subDir);

  // Book histograms
  booker.setCurrentFolder(subDir + "/HitProperty");
  clusterSize = booker.book1D("ClusterSize", "Cluster size;Cluster size", 11, -0.5, 10.5);
  clusterSizeBarrel = booker.book1D("ClusterSizeBarrel", "Cluster size in Barrel;Cluster size", 11, -0.5, 10.5);
  clusterSizeEndcap = booker.book1D("ClusterSizeEndcap", "Cluster size in Endcap;Cluster size", 11, -0.5, 10.5);

  avgClusterSize = booker.book1D("AverageClusterSize", "Average cluster size;Average clsuter size", 11, -0.5, 10.5);
  avgClusterSizeBarrel =
      booker.book1D("AverageClusterSizeBarrel", "Average cluster size in Barrel;Average clsuter size", 11, -0.5, 10.5);
  avgClusterSizeEndcap =
      booker.book1D("AverageClusterSizeEndcap", "Average cluster size in Endcap;Average clsuter size", 11, -0.5, 10.5);

  nRecHitBarrel =
      booker.book1D("NRecHitBarrel", "Number of RPC recHits per event in Barrel;Number of RPC hits", 25, 0, 25);
  nRecHitEndcap =
      booker.book1D("NRecHitEndcap", "Number of RPC recHits per event in Endcap;Number of RPC hits", 25, 0, 25);

  nRefHitBarrel =
      booker.book1D("NRefHitBarrel", "Number of reference hits per event in Barrel;Number of RPC hits", 25, 0, 25);
  nRefHitEndcap =
      booker.book1D("NRefHitEndcap", "Number of reference hits per event in Endcap;Number of RPC hits", 25, 0, 25);

  nMatchHitBarrel = booker.book1D(
      "nMatchBarrel", "Number of matched reference hits per event in Barrel;Number of RPC hits", 25, 0, 25);
  nMatchHitEndcap = booker.book1D(
      "nMatchEndcap", "Number of matched reference hits per event in Endcap;Number of RPC hits", 25, 0, 25);

  timeBarrel = booker.book1D("RecHitTimeBarrel", "RecHit time in Barrel;Time (ns)", 100, -12.5, 12.5);
  timeEndcap = booker.book1D("RecHitTimeEndcap", "RecHit time in Endcap;Time (ns)", 100, -12.5, 12.5);
  timeIRPC = booker.book1D("RecHitTimeIRPC", "RecHit time of iRPC;Time (ns)", 100, -12.5, 12.5);
  timeCRPC = booker.book1D("RecHitTimeCPRC", "RecHit time in current RPC;Time (ns)", 100, -12.5, 12.5);

  clusterSize->getTH1()->SetMinimum(0);
  clusterSizeBarrel->getTH1()->SetMinimum(0);
  clusterSizeEndcap->getTH1()->SetMinimum(0);

  avgClusterSize->getTH1()->SetMinimum(0);
  avgClusterSizeBarrel->getTH1()->SetMinimum(0);
  avgClusterSizeEndcap->getTH1()->SetMinimum(0);

  nRecHitBarrel->getTH1()->SetMinimum(0);
  nRecHitEndcap->getTH1()->SetMinimum(0);

  nRefHitBarrel->getTH1()->SetMinimum(0);
  nRefHitEndcap->getTH1()->SetMinimum(0);

  nMatchHitBarrel->getTH1()->SetMinimum(0);
  nMatchHitEndcap->getTH1()->SetMinimum(0);

  // Occupancy 1D
  booker.setCurrentFolder(subDir + "/Occupancy");
  refHitOccupancyBarrel_wheel = booker.book1D("RefHitOccupancyBarrel_wheel", "Reference Hit occupancy", 5, -2.5, 2.5);
  refHitOccupancyEndcap_disk = booker.book1D("RefHitOccupancyEndcap_disk", "Reference Hit occupancy", 9, -4.5, 4.5);
  refHitOccupancyBarrel_station =
      booker.book1D("RefHitOccupancyBarrel_station", "Reference Hit occupancy", 4, 0.5, 4.5);

  recHitOccupancyBarrel_wheel = booker.book1D("RecHitOccupancyBarrel_wheel", "RecHit occupancy", 5, -2.5, 2.5);
  recHitOccupancyEndcap_disk = booker.book1D("RecHitOccupancyEndcap_disk", "RecHit occupancy", 9, -4.5, 4.5);
  recHitOccupancyBarrel_station = booker.book1D("RecHitOccupancyBarrel_station", "RecHit occupancy", 4, 0.5, 4.5);

  matchOccupancyBarrel_wheel = booker.book1D("MatchOccupancyBarrel_wheel", "Matched hit occupancy", 5, -2.5, 2.5);
  matchOccupancyEndcap_disk = booker.book1D("MatchOccupancyEndcap_disk", "Matched hit occupancy", 9, -4.5, 4.5);
  matchOccupancyBarrel_station = booker.book1D("MatchOccupancyBarrel_station", "Matched hit occupancy", 4, 0.5, 4.5);

  refHitOccupancyBarrel_wheel->getTH1()->SetMinimum(0);
  refHitOccupancyEndcap_disk->getTH1()->SetMinimum(0);
  refHitOccupancyBarrel_station->getTH1()->SetMinimum(0);

  recHitOccupancyBarrel_wheel->getTH1()->SetMinimum(0);
  recHitOccupancyEndcap_disk->getTH1()->SetMinimum(0);
  recHitOccupancyBarrel_station->getTH1()->SetMinimum(0);

  matchOccupancyBarrel_wheel->getTH1()->SetMinimum(0);
  matchOccupancyEndcap_disk->getTH1()->SetMinimum(0);
  matchOccupancyBarrel_station->getTH1()->SetMinimum(0);

  // Occupancy 2D
  refHitOccupancyBarrel_wheel_station =
      booker.book2D("RefHitOccupancyBarrel_wheel_station", "Reference hit occupancy", 5, -2.5, 2.5, 4, 0.5, 4.5);
  refHitOccupancyEndcap_disk_ring =
      booker.book2D("RefHitOccupancyEndcap_disk_ring", "Reference hit occupancy", 9, -4.5, 4.5, 4, 0.5, 4.5);

  recHitOccupancyBarrel_wheel_station =
      booker.book2D("RecHitOccupancyBarrel_wheel_station", "RecHit occupancy", 5, -2.5, 2.5, 4, 0.5, 4.5);
  recHitOccupancyEndcap_disk_ring =
      booker.book2D("RecHitOccupancyEndcap_disk_ring", "RecHit occupancy", 9, -4.5, 4.5, 4, 0.5, 4.5);

  matchOccupancyBarrel_wheel_station =
      booker.book2D("MatchOccupancyBarrel_wheel_station", "Matched hit occupancy", 5, -2.5, 2.5, 4, 0.5, 4.5);
  matchOccupancyEndcap_disk_ring =
      booker.book2D("MatchOccupancyEndcap_disk_ring", "Matched hit occupancy", 9, -4.5, 4.5, 4, 0.5, 4.5);

  refHitOccupancyBarrel_wheel_station->getTH2F()->SetMinimum(0);
  refHitOccupancyEndcap_disk_ring->getTH2F()->SetMinimum(0);

  recHitOccupancyBarrel_wheel_station->getTH2F()->SetMinimum(0);
  recHitOccupancyEndcap_disk_ring->getTH2F()->SetMinimum(0);

  matchOccupancyBarrel_wheel_station->getTH2F()->SetMinimum(0);
  matchOccupancyEndcap_disk_ring->getTH2F()->SetMinimum(0);

  // Residuals
  booker.setCurrentFolder(subDir + "/Residual");
  resBarrel = booker.book1D("ResBarrel", "Global Residuals for Barrel;Residual [cm]", 100, -8, 8);
  resEndcap = booker.book1D("ResEndcap", "Global Residuals for Endcap;Residual [cm]", 100, -8, 8);

  resBarrel_W0 = booker.book1D("ResBarrel_W0", "Residual Wheel   0;Residual [cm]", 80, -8, 8);
  resBarrel_W1 = booker.book1D("ResBarrel_W1", "Residual Wheel +-1;Residual [cm]", 80, -8, 8);
  resBarrel_W2 = booker.book1D("ResBarrel_W2", "Residual Wheel +-2;Residual [cm]", 80, -8, 8);

  resBarrel_S1L1 = booker.book1D("ResBarrel_S1L1", "Residual Wheel Station1 Layer1;Residual [cm]", 80, -8, 8);
  resBarrel_S1L2 = booker.book1D("ResBarrel_S1L2", "Residual Wheel Station1 Layer2;Residual [cm]", 80, -8, 8);
  resBarrel_S2L1 = booker.book1D("ResBarrel_S2L1", "Residual Wheel Station2 Layer1;Residual [cm]", 80, -8, 8);
  resBarrel_S2L2 = booker.book1D("ResBarrel_S2L2", "Residual Wheel Station2 Layer2;Residual [cm]", 80, -8, 8);
  resBarrel_S3 = booker.book1D("ResBarrel_S3", "Residual Wheel Station3;Residual [cm]", 80, -8, 8);
  resBarrel_S4 = booker.book1D("ResBarrel_S4", "Residual Wheel Station4;Residual [cm]", 80, -8, 8);

  resEndcap_D1 = booker.book1D("ResEndcap_D1", "Residual Disk +-1;Residual [cm]", 80, -8, 8);
  resEndcap_D2 = booker.book1D("ResEndcap_D2", "Residual Disk +-2;Residual [cm]", 80, -8, 8);
  resEndcap_D3 = booker.book1D("ResEndcap_D3", "Residual Disk +-3;Residual [cm]", 80, -8, 8);
  resEndcap_D4 = booker.book1D("ResEndcap_D4", "Residual Disk +-4;Residual [cm]", 80, -8, 8);

  resEndcap_R1 = booker.book1D("ResEndcap_R1", "Residual Ring 1;Residual [cm]", 80, -8, 8);
  resEndcap_R2 = booker.book1D("ResEndcap_R2", "Residual Ring 2;Residual [cm]", 80, -8, 8);
  resEndcap_R3 = booker.book1D("ResEndcap_R3", "Residual Ring 3;Residual [cm]", 80, -8, 8);

  resBarrel->getTH1()->SetMinimum(0);
  resEndcap->getTH1()->SetMinimum(0);

  resBarrel_W0->getTH1()->SetMinimum(0);
  resBarrel_W0->getTH1()->SetMinimum(0);
  resBarrel_W0->getTH1()->SetMinimum(0);

  resBarrel_S1L1->getTH1()->SetMinimum(0);
  resBarrel_S1L2->getTH1()->SetMinimum(0);
  resBarrel_S2L1->getTH1()->SetMinimum(0);
  resBarrel_S2L2->getTH1()->SetMinimum(0);
  resBarrel_S3->getTH1()->SetMinimum(0);
  resBarrel_S4->getTH1()->SetMinimum(0);

  resEndcap_D1->getTH1()->SetMinimum(0);
  resEndcap_D2->getTH1()->SetMinimum(0);
  resEndcap_D3->getTH1()->SetMinimum(0);
  resEndcap_D4->getTH1()->SetMinimum(0);

  resEndcap_R1->getTH1()->SetMinimum(0);
  resEndcap_R2->getTH1()->SetMinimum(0);
  resEndcap_R3->getTH1()->SetMinimum(0);

  // Pulls
  pullBarrel = booker.book1D("PullBarrel", "Global Pull for Barrel;Pull", 100, -3, 3);
  pullEndcap = booker.book1D("PullEndcap", "Global Pull for Endcap;Pull", 100, -3, 3);

  pullBarrel_W0 = booker.book1D("ResBarrel_W0", "Residual Wheel   0;Residual [cm]", 80, -8, 8);
  pullBarrel_W1 = booker.book1D("ResBarrel_W1", "Residual Wheel +-1;Residual [cm]", 80, -8, 8);
  pullBarrel_W2 = booker.book1D("ResBarrel_W2", "Residual Wheel +-2;Residual [cm]", 80, -8, 8);

  pullBarrel_S1L1 = booker.book1D("ResBarrel_S1L1", "Residual Wheel Station1 Layer1;Residual [cm]", 80, -8, 8);
  pullBarrel_S1L2 = booker.book1D("ResBarrel_S1L2", "Residual Wheel Station1 Layer2;Residual [cm]", 80, -8, 8);
  pullBarrel_S2L1 = booker.book1D("ResBarrel_S2L1", "Residual Wheel Station2 Layer1;Residual [cm]", 80, -8, 8);
  pullBarrel_S2L2 = booker.book1D("ResBarrel_S2L2", "Residual Wheel Station2 Layer2;Residual [cm]", 80, -8, 8);
  pullBarrel_S3 = booker.book1D("ResBarrel_S3", "Residual Wheel Station3;Residual [cm]", 80, -8, 8);
  pullBarrel_S4 = booker.book1D("ResBarrel_S4", "Residual Wheel Station4;Residual [cm]", 80, -8, 8);

  pullEndcap_D1 = booker.book1D("ResEndcap_D1", "Residual Disk +-1;Residual [cm]", 80, -8, 8);
  pullEndcap_D2 = booker.book1D("ResEndcap_D2", "Residual Disk +-2;Residual [cm]", 80, -8, 8);
  pullEndcap_D3 = booker.book1D("ResEndcap_D3", "Residual Disk +-3;Residual [cm]", 80, -8, 8);
  pullEndcap_D4 = booker.book1D("ResEndcap_D4", "Residual Disk +-4;Residual [cm]", 80, -8, 8);

  pullEndcap_R1 = booker.book1D("ResEndcap_R1", "Residual Ring 1;Residual [cm]", 80, -8, 8);
  pullEndcap_R2 = booker.book1D("ResEndcap_R2", "Residual Ring 2;Residual [cm]", 80, -8, 8);
  pullEndcap_R3 = booker.book1D("ResEndcap_R3", "Residual Ring 3;Residual [cm]", 80, -8, 8);

  pullBarrel->getTH1()->SetMinimum(0);
  pullEndcap->getTH1()->SetMinimum(0);

  pullBarrel_W0->getTH1()->SetMinimum(0);
  pullBarrel_W0->getTH1()->SetMinimum(0);
  pullBarrel_W0->getTH1()->SetMinimum(0);

  pullBarrel_S1L1->getTH1()->SetMinimum(0);
  pullBarrel_S1L2->getTH1()->SetMinimum(0);
  pullBarrel_S2L1->getTH1()->SetMinimum(0);
  pullBarrel_S2L2->getTH1()->SetMinimum(0);
  pullBarrel_S3->getTH1()->SetMinimum(0);
  pullBarrel_S4->getTH1()->SetMinimum(0);

  pullEndcap_D1->getTH1()->SetMinimum(0);
  pullEndcap_D2->getTH1()->SetMinimum(0);
  pullEndcap_D3->getTH1()->SetMinimum(0);
  pullEndcap_D4->getTH1()->SetMinimum(0);

  pullEndcap_R1->getTH1()->SetMinimum(0);
  pullEndcap_R2->getTH1()->SetMinimum(0);
  pullEndcap_R3->getTH1()->SetMinimum(0);

  // Set plot options
  refHitOccupancyBarrel_wheel_station->getTH2F()->SetOption("COLZ");
  refHitOccupancyEndcap_disk_ring->getTH2F()->SetOption("COLZ");
  recHitOccupancyBarrel_wheel_station->getTH2F()->SetOption("COLZ");
  recHitOccupancyEndcap_disk_ring->getTH2F()->SetOption("COLZ");
  matchOccupancyBarrel_wheel_station->getTH2F()->SetOption("COLZ");
  matchOccupancyEndcap_disk_ring->getTH2F()->SetOption("COLZ");

  refHitOccupancyBarrel_wheel_station->getTH2F()->SetContour(10);
  refHitOccupancyEndcap_disk_ring->getTH2F()->SetContour(10);
  recHitOccupancyBarrel_wheel_station->getTH2F()->SetContour(10);
  recHitOccupancyEndcap_disk_ring->getTH2F()->SetContour(10);
  matchOccupancyBarrel_wheel_station->getTH2F()->SetContour(10);
  matchOccupancyEndcap_disk_ring->getTH2F()->SetContour(10);

  refHitOccupancyBarrel_wheel_station->getTH2F()->SetStats(false);
  refHitOccupancyEndcap_disk_ring->getTH2F()->SetStats(false);
  recHitOccupancyBarrel_wheel_station->getTH2F()->SetStats(false);
  recHitOccupancyEndcap_disk_ring->getTH2F()->SetStats(false);
  matchOccupancyBarrel_wheel_station->getTH2F()->SetStats(false);
  matchOccupancyEndcap_disk_ring->getTH2F()->SetStats(false);

  // Set bin labels
  for (int i = 1; i <= 5; ++i) {
    TString binLabel = Form("Wheel %d", i - 3);

    refHitOccupancyBarrel_wheel->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);
    recHitOccupancyBarrel_wheel->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);
    matchOccupancyBarrel_wheel->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);

    refHitOccupancyBarrel_wheel_station->getTH2F()->GetXaxis()->SetBinLabel(i, binLabel);
    recHitOccupancyBarrel_wheel_station->getTH2F()->GetXaxis()->SetBinLabel(i, binLabel);
    matchOccupancyBarrel_wheel_station->getTH2F()->GetXaxis()->SetBinLabel(i, binLabel);
  }

  for (int i = 1; i <= 9; ++i) {
    TString binLabel = Form("Disk %d", i - 5);

    refHitOccupancyEndcap_disk->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);
    recHitOccupancyEndcap_disk->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);
    matchOccupancyEndcap_disk->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);

    refHitOccupancyEndcap_disk_ring->getTH2F()->GetXaxis()->SetBinLabel(i, binLabel);
    recHitOccupancyEndcap_disk_ring->getTH2F()->GetXaxis()->SetBinLabel(i, binLabel);
    matchOccupancyEndcap_disk_ring->getTH2F()->GetXaxis()->SetBinLabel(i, binLabel);
  }

  for (int i = 1; i <= 4; ++i) {
    TString binLabel = Form("Station %d", i);

    refHitOccupancyBarrel_station->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);
    recHitOccupancyBarrel_station->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);
    matchOccupancyBarrel_station->getTH1F()->GetXaxis()->SetBinLabel(i, binLabel);

    refHitOccupancyBarrel_wheel_station->getTH2F()->GetYaxis()->SetBinLabel(i, binLabel);
    recHitOccupancyBarrel_wheel_station->getTH2F()->GetYaxis()->SetBinLabel(i, binLabel);
    matchOccupancyBarrel_wheel_station->getTH2F()->GetYaxis()->SetBinLabel(i, binLabel);
  }

  for (int i = 1; i <= 4; ++i) {
    TString binLabel = Form("Ring %d", i);

    refHitOccupancyEndcap_disk_ring->getTH2F()->GetYaxis()->SetBinLabel(i, binLabel);
    recHitOccupancyEndcap_disk_ring->getTH2F()->GetYaxis()->SetBinLabel(i, binLabel);
    matchOccupancyEndcap_disk_ring->getTH2F()->GetYaxis()->SetBinLabel(i, binLabel);
  }

  booked_ = true;

  booker.setCurrentFolder(pwd);
}
