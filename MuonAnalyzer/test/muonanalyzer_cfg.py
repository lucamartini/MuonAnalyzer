import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# import of standard configurations
# process.load('Configuration.StandardSequences.Services_cff')
# process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
# process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
# process.load('Configuration.StandardSequences.RawToDigi_cff')
# process.load('Configuration.StandardSequences.L1Reco_cff')
# process.load('Configuration.StandardSequences.Reconstruction_cff')
# process.load('CommonTools.ParticleFlow.EITopPAG_cff')
# process.load('Configuration.StandardSequences.Validation_cff')
# process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
# process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_74_V9::All', '')


process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(
                                '/store/relval/CMSSW_7_4_6/RelValBuJpsiK_13/GEN-SIM-RECO/MCRUN2_74_V9-v2/00000/50D0E1C6-2B1A-E511-A9B4-0025905A608C.root',
                                '/store/relval/CMSSW_7_4_6/RelValBuJpsiK_13/GEN-SIM-RECO/MCRUN2_74_V9-v2/00000/D6EA9003-2C1A-E511-BCE0-0025905A60C6.root',
                            )
)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(100))


process.demo = cms.EDAnalyzer('MuonAnalyzer',
                              OutputFileName = cms.string("muontree.root"),
                              doReco = cms.bool(True),
                              doMC = cms.bool(True),
                              GenParticleCollection = cms.InputTag("genParticles"),
                              doTrigger = cms.bool(True),
                              doRAWTrigger = cms.bool(False),
                              pdgId = cms.vint32(13, 433),
                              HLTString = cms.string("HLT"),
                              HLTPaths = cms.vstring(
        "HLT_Dimuon0_Jpsi_Muon_v",
        "HLT_Dimuon0er16_Jpsi_NoOS_NoVertexing_v",
        "HLT_Dimuon0er16_Jpsi_NoVertexing_v",
        "HLT_Dimuon10_Jpsi_Barrel_v",
        "HLT_Dimuon13_PsiPrime_v",
        "HLT_Dimuon16_Jpsi_v",
        "HLT_Dimuon20_Jpsi_v",
        "HLT_Dimuon6_Jpsi_NoVertexing_v",
        "HLT_Dimuon8_PsiPrime_Barrel_v",
        "HLT_DoubleMu4_3_Bs_v",
        "HLT_DoubleMu4_3_Jpsi_Displaced_v",
        "HLT_DoubleMu4_JpsiTrk_Displaced_v",
        "HLT_DoubleMu4_PsiPrimeTrk_Displaced_v",
        "HLT_Mu7p5_L2Mu2_Jpsi_v",
        "HLT_Mu7p5_Track2_Jpsi_v",
        "HLT_Mu7p5_Track3p5_Jpsi_v",
        "HLT_Mu7p5_Track7_Jpsi_v",
        "HLT_QuadMuon0_Dimuon0_Jpsi_v",
        "HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v",
        "HLT_Dimuon0_Phi_Barrel_v",
        "HLT_Dimuon0_Upsilon_Muon_v",
        "HLT_Dimuon13_Upsilon_v",
        "HLT_Dimuon8_Upsilon_Barrel_v",
        "HLT_Mu16_TkMu0_dEta18_Onia_v",
        "HLT_Mu16_TkMu0_dEta18_Phi_v",
        "HLT_Mu25_TkMu0_dEta18_Onia_v",
        "HLT_Mu7p5_L2Mu2_Upsilon_v",
        "HLT_Mu7p5_Track2_Upsilon_v",
        "HLT_Mu7p5_Track3p5_Upsilon_v",
        "HLT_Mu7p5_Track7_Upsilon_v",
        "HLT_QuadMuon0_Dimuon0_Upsilon_v",


        ),
                              OnlineBeamSpot = cms.InputTag("onlineBeamSpot"),
                              OfflineBeamSpot = cms.InputTag("offlineBeamSpot"),
                              PrimaryVertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.p = cms.Path(process.demo)

# customisation of the process.


# End of customisation functions
