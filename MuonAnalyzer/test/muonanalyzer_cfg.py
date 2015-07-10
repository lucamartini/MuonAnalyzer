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
                              HLTPaths = cms.vstring("HLT_DoubleMu4_3_Jpsi_Displaced_v"),
                              OnlineBeamSpot = cms.InputTag("onlineBeamSpot"),
                              OfflineBeamSpot = cms.InputTag("offlineBeamSpot"),
                              PrimaryVertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.p = cms.Path(process.demo)

# customisation of the process.


# End of customisation functions
