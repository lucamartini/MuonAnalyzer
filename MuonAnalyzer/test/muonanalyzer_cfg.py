import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(20) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_100_1_lhU.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_101_1_PgD.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_10_1_aGo.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_11_1_XUW.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_12_1_fjY.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_13_1_xSj.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_14_1_aWs.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_15_1_2fk.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_16_1_SmM.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_17_1_q8V.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_18_1_Kmi.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_19_1_eex.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_1_1_qZa.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_20_1_WFS.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_21_1_JrL.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_22_1_cg4.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_23_1_S5r.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_24_1_CVg.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_25_1_kfe.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_26_1_Xur.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_27_1_tzd.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_28_1_ZIC.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_29_1_4n0.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_2_1_jAf.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_30_1_PVm.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_31_1_HJU.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_32_1_WRf.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_33_1_tME.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_34_1_SXn.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_35_1_kC3.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_36_1_97h.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_37_1_7tV.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_38_1_nVL.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_39_1_o4i.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_3_1_9g3.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_40_1_9dj.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_41_1_b2N.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_42_1_zyB.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_43_1_zE9.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_44_1_TGy.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_45_1_GyY.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_46_1_OwU.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_47_1_QdP.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_48_1_9g5.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_49_1_PX5.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_4_1_nJE.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_50_1_Gwx.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_51_1_H4a.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_52_1_nBh.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_53_1_9kD.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_54_1_OlC.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_55_1_npE.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_56_1_pG9.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_57_1_hBS.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_58_1_82S.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_59_1_owc.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_5_1_2o0.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_60_1_cu8.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_61_1_tu4.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_62_1_rNp.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_63_1_G5n.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_64_1_MuO.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_66_1_iQz.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_67_1_B2Z.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_68_1_YTZ.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_69_1_sdH.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_70_1_7aR.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_71_1_Cte.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_72_1_umV.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_73_1_qQr.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_74_1_Leg.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_75_1_cl1.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_76_1_9oz.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_77_1_Nj4.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_78_1_4Lz.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_79_1_GVr.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_7_1_rm3.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_80_1_Bj3.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_81_1_WGS.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_82_1_sZY.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_83_1_IwZ.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_84_1_fKE.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_85_1_mie.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_86_1_Mrj.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_87_1_gsF.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_88_1_nNs.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_89_1_uwW.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_8_1_rUu.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_90_1_0YJ.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_91_1_VcL.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_92_1_OSV.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_93_1_TRp.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_94_1_emH.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_95_1_PWu.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_96_1_og8.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_97_1_Qoj.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_98_1_0HY.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_99_1_Reb.root',
'/store/user/lmartini/BsMM_5314_GEN-SIM_v4/BsMM_5314_RAW2DIGI_L1Reco_RECO_v4/e42f2f19740e5f379fd108eac1484ecb/BsMM_step3_5314_9_1_8MR.root',
    )
)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')

process.demo = cms.EDAnalyzer('MuonAnalyzer',
    OutputFileName = cms.string("muontree_5314_v2.root"),
    doMC = cms.bool(True),
    doReco = cms.bool(True),
    doTrigger = cms.bool(True),
    doRAWTrigger = cms.bool(False),
    pdgId = cms.vint32(13, 531),
    HLTString = cms.string("HLT"),
    HLTPaths = cms.vstring("HLT_DoubleMu4_Jpsi_Displaced_v", "HLT_DoubleMu4_JpsiTk_Displaced_v", "HLT_Dimuon0_Jpsi_v", "HLT_Dimuon0_Jpsi_NoVertexing_v", "HLT_Dimuon8_Jpsi_v", "HLT_Dimuon10_Jpsi_v", "HLT_Dimuon0_Jpsi_Muon_v", "HLT_Mu5_L2Mu3_Jpsi_v", "HLT_Mu5_Track2_Jpsi_v22", "HLT_Mu5_Track3p5_Jpsi_v", "HLT_Mu7_Track7_Jpsi_v")
)


process.p = cms.Path(process.demo)
