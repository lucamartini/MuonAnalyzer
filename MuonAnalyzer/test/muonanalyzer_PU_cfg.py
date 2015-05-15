import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# import of standard configurations
# process.load('Configuration.StandardSequences.Services_cff')
# process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mix_POISSON_average_cfi')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
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
process.GlobalTag = GlobalTag(process.GlobalTag, 'PH2_1K_FB_V2::All', '')


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/02BFD693-A1DE-E311-BFF0-00248C0BE014.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/106F1889-EADC-E311-8FB8-002590596486.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/12007767-57DC-E311-9311-0025905A48E4.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/128E0E90-84DE-E311-BB1F-002618943949.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/14CD8821-CEDE-E311-A961-0025905A60A8.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/16E3C529-D7DC-E311-BBAD-0025905A6134.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/1C69241A-29DD-E311-81C1-0025905A609A.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/1EF61585-5CDC-E311-ADA1-003048FFCBA4.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/1EF9737E-C0DC-E311-ADF5-0025905964B4.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/2879C3F7-F4DC-E311-8DA5-002618943967.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/40D36320-82DD-E311-9D73-003048679076.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/4464C91F-75DC-E311-9741-003048FFD7BE.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/487B54AF-07DD-E311-82EE-002590593902.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/50CC659A-57DD-E311-901D-0025905AA9F0.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/5645C5A4-A5DD-E311-BA37-0025905A6066.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/5AB67BF8-63DC-E311-A8EB-0025905A60A6.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/62766C29-66DC-E311-93F1-0025905A48D8.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/6485C729-C7DC-E311-ACFC-0025905A60D2.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/64AF398E-B7DC-E311-8F52-002618943821.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/68ADA9DA-70DD-E311-B7FA-0025905A6092.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/6AC06C4B-60DC-E311-BC40-003048FFD740.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/724B2604-6ADC-E311-AC05-0025905A48D8.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/78F833A0-11DD-E311-B752-003048FFD76E.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/847B3FA7-9DDC-E311-B8FA-002618943821.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/86CA706D-1DDD-E311-BAF4-00261894387E.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/980960DA-93DC-E311-B79F-0025905938AA.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/A0A729D9-BADD-E311-BCEC-002354EF3BDB.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/AC306989-12DE-E311-B198-002618943821.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/B03E2A29-77DE-E311-B7AE-002354EF3BE0.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/B21C5073-ABDC-E311-A413-002618943975.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/B4D7307B-83DC-E311-9F9D-0025905A6090.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/C24A02BF-7CDE-E311-8F39-002354EF3BDF.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/C45C534C-4ADD-E311-9AE8-0025905A612C.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/CADAB0CC-FCDC-E311-9FA4-0026189438FE.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/DE136C97-7FDC-E311-A291-0025905A60EE.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/E6062998-98DD-E311-BE99-002590596498.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/EE476C5F-3BDD-E311-9DC9-0025905A60A6.root",
"/store/mc/Muon2023Upg14DR/Bs2MuMu_TuneZ2star_2023_14TeV_Pythia6/AODSIM/PU140bx25_PH2_1K_FB_V2-v1/00000/FAE4E17A-B1DC-E311-BC69-00261894395F.root",
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


process.demo = cms.EDAnalyzer('MuonAnalyzer',
                              OutputFileName = cms.string("muontree_620SLHC7_PU.root"),
                              doReco = cms.bool(True),
                              doMC = cms.bool(True),
                              doTrigger = cms.bool(False),
                              doRAWTrigger = cms.bool(False),
                              pdgId = cms.vint32(13, 531),
                              HLTString = cms.string("HLT"),
                              HLTPaths = cms.vstring("HLT_DoubleMu4_Jpsi_Displaced_v", "HLT_Mu7_Track7_Jpsi_v"),
                              OnlineBeamSpot = cms.InputTag("onlineBeamSpot"),
                              OfflineBeamSpot = cms.InputTag("offlineBeamSpot"),
                              PrimaryVertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.p = cms.Path(process.demo)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023Muon

#call to customisation function cust_2023Muon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023Muon(process)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# End of customisation functions
