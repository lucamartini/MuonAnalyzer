import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12_DR53X/BsToMuMu_8TeV-pythia6-evtgen/AODSIM/PU_S10_START53_V7C-v1/10000/FEF1BB3C-EB51-E211-944A-00304867C0EA.root',
    )
)

process.demo = cms.EDAnalyzer('MuonAnalyzer',
    OutputFileName = cms.string("muontree.root"),
    IsMC = cms.bool(True),
    pdgId = cms.vint32(13, 13),
)


process.p = cms.Path(process.demo)
