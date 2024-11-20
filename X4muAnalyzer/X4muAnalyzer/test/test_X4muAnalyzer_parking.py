import FWCore.ParameterSet.Config as cms
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process('X4mu')
#process = cms.Process("X4muAnalyzer")
ouput_filename = 'X4mu_parking.root'

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_PromptAnalysis_v1', '') #2022 ABCDE PromptReco
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_v15', '') #2022 ABCDE ReReco
#process.GlobalTag = GlobalTag(process.GlobalTag, '124X_dataRun3_PromptAnalysis_v2', '') #2022 FG PromptReco
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_dataRun3_PromptAnalysis_v1', '') #2023 CD ReReco

process.MessageLogger.cerr.FwkSummary.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) # -1 
)

process.source = cms.Source("PoolSource",
        # replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024H/ParkingDoubleMuonLowMass1/MINIAOD/PromptReco-v1/000/385/836/00000/3c074159-71ac-450d-9b91-e2eb2c2a265e.root'
        )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ouput_filename),
)

process.options = cms.untracked.PSet(
    TryToContinue = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    Rethrow = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(8),
    numberOfThreads = cms.untracked.uint32(8),
    printDependencies = cms.untracked.bool(False),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False),
)

process.X4muMuonFilter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("slimmedMuons"),
    minNumber = cms.uint32(2),
)

process.X4muTrackFilter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("packedPFCandidates"),
    minNumber = cms.uint32(6),
)

process.X4muVertexFinder = cms.EDProducer("X4muPatSecondaryVertexProducer",
    tracks = cms.InputTag("packedPFCandidates"),
    patMuon = cms.InputTag("slimmedMuons"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    MesonPaiMass = cms.double(2.9841), # Eta_c mass
    MesonPaiMassErr = cms.double(0.0004), # Eta_c mass error
    MesonMuMass = cms.double(3.0969), # J/psi mass
    MesonMuMassErr = cms.double(0.00004), # J/psi mass error
    MesonMassWindow = cms.double(0.5), # 500 MeV
    ExMesonMass = cms.double(-999.0), # Not used
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'HLT_DoubleMu2_Jpsi_LowPt',                #  1=        1
    'HLT_Dimuon0_Jpsi3p5_Muon2',               #  2=        2
    'HLT_DoubleMu4_JpsiTrkTrk_Displaced',      #  3=        4
    'HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05'    #  4=        8
    ),
    doIso = cms.bool(True),
    MuIso = cms.double(0.02),
    PaiIso = cms.double(0.02),
    doTrigger = cms.bool(False),
    vProb = cms.double(0.01),
    selectionType = cms.int32(0), # Not used
    maxLoop = cms.uint32(1000000), # Max Loop
)

process.p = cms.Path(process.X4muMuonFilter*process.X4muTrackFilter*process.X4muVertexFinder)