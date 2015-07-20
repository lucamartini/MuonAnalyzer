import FWCore.ParameterSet.Config as cms

process = cms.Process("Muonanalyzer")

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


process.load("MuonAnalyzer.MuonAnalyzer.MuonAnalyzer_cfi")

process.p = cms.Path(process.muonanalyzer)


