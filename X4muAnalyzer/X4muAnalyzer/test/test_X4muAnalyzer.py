import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer
from Configuration.Eras.Era_Run3_2024_cff import Run3_2024
from Configuration.AlCa.GlobalTag import GlobalTag

process = cms.Process('NANO', Run3_2024)
#process = cms.Process("X4muAnalyzer")
ouput_filename = 'X4mu_scouting.root'

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('PhysicsTools.NanoAOD.nano_cff')
process.load('PhysicsTools.NanoAOD.common_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")

process.GlobalTag = GlobalTag(process.GlobalTag, '140X_dataRun3_Prompt_v4', '')

process.MessageLogger.cerr.FwkSummary.reportEvery = 10000
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) # -1 
)

process.source = cms.Source("PoolSource",
        # replace 'myfile.root' with the source file you want to use
        fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch/' + '/store/data/Run2024G/ScoutingPFRun3/HLTSCOUT/v1/000/385/312/00000/72e81dd1-3997-44ad-a166-dbc0bbb1d1f4.root'
        )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(ouput_filename),
    #fileOption = cms.string('RECREATE'),
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring("ProductNotFound","TooManyProducts","TooFewProducts"),
    TryToContinue = cms.untracked.vstring('ProductNotFound'),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(1),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False),
)

process.X4muConverter = cms.EDProducer("X4muScoutingToRecoMuonProducer",
    scoutingMuon = cms.InputTag("hltScoutingMuonPackerVtx"),
    scoutingMuonNoVtx = cms.InputTag("hltScoutingMuonPackerNoVtx"),
    scoutingTrack = cms.InputTag("hltScoutingTrackPacker"),
)

process.X4muFilter = cms.EDFilter('CandViewCountFilter',
    src       = cms.InputTag("X4muConverter", "recoMuons"),
    minNumber = cms.uint32(4),  # 4
)

process.X4muVertexFinder = cms.EDProducer("X4muSecondaryVertexProducer",
    recoMuon = cms.InputTag("X4muConverter", "recoMuons"),
    recoTrack = cms.InputTag("X4muConverter", "recoTracks"),
    MesonMassBig = cms.double(3.0969), # J/Psi
    MesonMassBigErr = cms.double(0.00004), # J/Psi
    MesonMassSmall = cms.double(3.0969), # J/Psi
    MesonMassSmallErr = cms.double(0.00004), # J/Psi
)

#MesonMassBig = cms.double(9.4604), # Upsilon
#MesonMassBigErr = cms.double(0.0001), # Upsilon

process.p = cms.Path(process.X4muConverter*process.X4muFilter*process.X4muVertexFinder)