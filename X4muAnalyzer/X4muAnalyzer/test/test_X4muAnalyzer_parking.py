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
#process.load('Configuration.StandardSequences.Skims_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v37', '')

process.MessageLogger.cerr.FwkSummary.reportEvery = 1000000000
process.MessageLogger.cerr.FwkReport.reportEvery = 1000000000

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) # -1 
)

process.source = cms.Source("PoolSource",
        # replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch/' + '/store/data/Run2016B/MuOnia/MINIAOD/ver1_HIPM_UL2016_MiniAODv2-v1/120000/D1125B17-272D-C945-8466-68CC548D4216.root'
        )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ouput_filename),
    #fileOption = cms.string('RECREATE'),
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    Rethrow = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    canDeleteEarly = cms.untracked.vstring(),
    fileMode = cms.untracked.string('FULLMERGE'),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False),
)

process.X4muFilter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("slimmedMuons"),
    minNumber = cms.uint32(4),  # 4
)

process.X4muVertexFinder = cms.EDProducer("X4muPatSecondaryVertexProducer",
    patMuon = cms.InputTag("slimmedMuons"),
    MesonMass = cms.double(9.4604), # Upsilon 9.4604
    MesonMassErr = cms.double(0.0001), # Upsilon 0.0001
    ExMesonMass = cms.double(3.0969), # J/Psi 3.0969
    TriggerResults  = cms.InputTag("TriggerResults", "", "HLT"),
    FilterNames     = cms.vstring(
    'HLT_Dimuon0_Upsilon_Muon',       #  1=        1
    "HLT_Dimuon0_Jpsi_Muon_v"         #  2=        2
    ),
    doIso = cms.bool(True),
)

process.p = cms.Path(process.X4muFilter*process.X4muVertexFinder)