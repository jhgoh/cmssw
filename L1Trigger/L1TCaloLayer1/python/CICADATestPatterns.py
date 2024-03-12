import FWCore.ParameterSet.Config as cms

standardCICADATestPatterns = cms.VPSet(
    cms.PSet(
        iPhi_1 = cms.vuint32(0, 0, 1, 0, 1, 0, 2, 3, 0, 0, 0, 3, 6, 0, ),
        iPhi_2 = cms.vuint32(2, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 2, ),
        iPhi_3 = cms.vuint32(0, 0, 0, 0, 0, 1, 0, 0, 5, 0, 0, 0, 0, 1, ),
        iPhi_4 = cms.vuint32(0, 1, 0, 0, 0, 1, 0, 0, 31, 1, 8, 7, 2, 8, ),
        iPhi_5 = cms.vuint32(1, 0, 1, 0, 0, 1, 0, 1, 2, 4, 0, 0, 0, 0, ),
        iPhi_6 = cms.vuint32(0, 0, 0, 0, 4, 0, 0, 1, 0, 0, 0, 0, 0, 6, ),
        iPhi_7 = cms.vuint32(0, 3, 1, 2, 1, 5, 1, 0, 0, 0, 0, 0, 1, 1, ),
        iPhi_8 = cms.vuint32(0, 0, 3, 2, 0, 2, 3, 3, 8, 10, 1, 2, 0, 27, ),
        iPhi_9 = cms.vuint32(6, 0, 0, 2, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, ),
        iPhi_10 = cms.vuint32(0, 0, 1, 0, 12, 2, 0, 0, 0, 1, 0, 1, 0, 2, ),
        iPhi_11 = cms.vuint32(5, 0, 0, 1, 0, 0, 1, 4, 2, 0, 15, 0, 0, 212, ),
        iPhi_12 = cms.vuint32(4, 0, 2, 0, 2, 1, 1, 4, 1, 0, 2, 3, 0, 0, ),
        iPhi_13 = cms.vuint32(0, 4, 1, 2, 182, 0, 2, 2, 0, 0, 0, 1, 1, 0, ),
        iPhi_14 = cms.vuint32(0, 10, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 2, ),
        iPhi_15 = cms.vuint32(6, 1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 12, ),
        iPhi_16 = cms.vuint32(0, 0, 0, 1, 0, 1, 0, 0, 3, 1, 0, 0, 0, 1, ),
        iPhi_17 = cms.vuint32(0, 0, 0, 0, 0, 2, 0, 4, 2, 0, 3, 0, 0, 2, ),
        iPhi_18 = cms.vuint32(2, 0, 0, 0, 0, 1, 0, 4, 0, 2, 4, 5, 0, 0, ),                                                
    ),
    cms.PSet(
        iPhi_1 = cms.vuint32(3, 5, 6, 2, 1, 0, 9, 0, 1, 1, 2, 1, 1, 5, ),
        iPhi_2 = cms.vuint32(4, 2, 1, 0, 5, 0, 0, 2, 4, 11, 10, 1, 1, 12, ),
        iPhi_3 = cms.vuint32(5, 0, 0, 2, 1, 2, 1, 1, 19, 20, 237, 0, 2, 2, ),
        iPhi_4 = cms.vuint32(5, 1, 0, 3, 2, 1, 2, 3, 3, 1, 2, 1, 1, 7, ),
        iPhi_5 = cms.vuint32(1, 1, 1, 2, 0, 0, 0, 3, 5, 2, 1, 1, 3, 14, ),
        iPhi_6 = cms.vuint32(4, 0, 2, 2, 0, 0, 0, 2, 1, 3, 3, 1, 0, 3, ),
        iPhi_7 = cms.vuint32(1, 4, 62, 6, 0, 1, 10, 2, 2, 5, 1, 1, 0, 7, ),
        iPhi_8 = cms.vuint32(13, 1, 0, 2, 1, 5, 1, 3, 1, 0, 1, 0, 4, 2, ),
        iPhi_9 = cms.vuint32(4, 1, 2, 1, 6, 2, 6, 0, 2, 2, 1, 0, 0, 6, ),
        iPhi_10 = cms.vuint32(10, 0, 2, 0, 3, 0, 1, 2, 12, 0, 20, 4, 0, 7, ),
        iPhi_11 = cms.vuint32(16, 2, 4, 1, 0, 2, 3, 15, 4, 1, 0, 6, 5, 5, ),
        iPhi_12 = cms.vuint32(3, 0, 1, 0, 1, 1, 4, 2, 9, 115, 38, 2, 3, 1, ),
        iPhi_13 = cms.vuint32(10, 3, 10, 15, 2, 0, 8, 8, 0, 2, 2, 0, 1, 8, ),
        iPhi_14 = cms.vuint32(4, 0, 0, 0, 1, 4, 0, 1, 1, 1, 1, 1, 0, 2, ),
        iPhi_15 = cms.vuint32(11, 1, 1, 2, 1, 3, 5, 4, 4, 2, 0, 1, 0, 13, ),
        iPhi_16 = cms.vuint32(6, 1, 1, 1, 0, 1, 3, 2, 1, 10, 3, 0, 0, 15, ),
        iPhi_17 = cms.vuint32(4, 0, 0, 1, 2, 1, 1, 2, 0, 1, 0, 1, 0, 3, ),
        iPhi_18 = cms.vuint32(5, 0, 0, 0, 4, 1, 0, 2, 5, 31, 0, 1, 1, 5, ),
    ),
    cms.PSet(
        iPhi_1 = cms.vuint32(4, 2, 2, 0, 0, 0, 4, 6, 1, 0, 0, 2, 2, 7, ),
        iPhi_2 = cms.vuint32(2, 2, 0, 1, 1, 1, 0, 0, 1, 2, 2, 1, 0, 0, ),
        iPhi_3 = cms.vuint32(0, 0, 0, 1, 52, 0, 3, 2, 7, 2, 0, 0, 1, 4, ),
        iPhi_4 = cms.vuint32(4, 0, 0, 0, 51, 6, 53, 4, 1, 0, 0, 0, 0, 0, ),
        iPhi_5 = cms.vuint32(10, 0, 0, 0, 1, 1, 4, 1, 0, 0, 0, 0, 0, 8, ),
        iPhi_6 = cms.vuint32(2, 0, 2, 0, 1, 5, 1, 3, 4, 0, 1, 0, 1, 14, ),
        iPhi_7 = cms.vuint32(1, 0, 1, 1, 0, 0, 8, 9, 2, 3, 0, 1, 0, 3, ),
        iPhi_8 = cms.vuint32(4, 0, 23, 62, 31, 0, 5, 3, 3, 1, 0, 0, 0, 4, ),
        iPhi_9 = cms.vuint32(100, 3, 10, 5, 0, 2, 0, 2, 1, 2, 0, 0, 0, 0, ),
        iPhi_10 = cms.vuint32(27, 2, 0, 0, 0, 2, 3, 1, 3, 0, 0, 2, 0, 0, ),
        iPhi_11 = cms.vuint32(8, 2, 3, 5, 5, 1, 1, 0, 4, 2, 2, 0, 0, 5, ),
        iPhi_12 = cms.vuint32(6, 6, 1, 0, 0, 2, 0, 3, 1, 3, 2, 1, 0, 2, ),
        iPhi_13 = cms.vuint32(0, 2, 2, 1, 0, 0, 7, 6, 0, 0, 0, 0, 1, 352, ),
        iPhi_14 = cms.vuint32(8, 0, 0, 1, 1, 1, 2, 2, 1, 4, 0, 0, 0, 2, ),
        iPhi_15 = cms.vuint32(3, 0, 0, 0, 1, 3, 3, 3, 0, 1, 0, 0, 0, 2, ),
        iPhi_16 = cms.vuint32(3, 166, 0, 4, 0, 2, 3, 1, 1, 1, 0, 0, 0, 6, ),
        iPhi_17 = cms.vuint32(2, 2, 1, 0, 0, 0, 0, 2, 5, 0, 0, 0, 0, 2, ),
        iPhi_18 = cms.vuint32(6, 3, 0, 2, 0, 4, 7, 1, 4, 4, 0, 0, 1, 2, ),
    ),
)
