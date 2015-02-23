import FWCore.ParameterSet.Config as cms

process = cms.Process("IIHEAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
'/store/mc/TP2023SHCALDR/ZprimeSSMToEE_M-1000_TuneZ2star_14TeV-pythia6/GEN-SIM-RECO/SHCALJan23_PU140BX25_PH2_1K_FB_V6-v1/20000/06E0FFCC-C1A4-E411-8821-008CFA000BB8.root'
])

# Global tags:
# PHYS14 25ns: PHYS14_25_V1

globalTag = 'PH2_1K_FB_V6::All'
process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag

process.out = cms.OutputModule("PoolOutputModule",
    ##process.FEVTSIMEventContent,
    fileName = cms.untracked.string('outfile.root')
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100000)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('outfile.root')
)

process.load("UserCode.IIHETree.IIHETree_cfi")
process.IIHEAnalysis.globalTag = cms.string(globalTag)

# Set pt or mass thresholds for the truth module here
# Setting thresholds reduces the size of the output files significantly 
process.IIHEAnalysis.MCTruth_ptThreshold = cms.untracked.double(10.0)
process.IIHEAnalysis.MCTruth_mThreshold  = cms.untracked.double(20.0)

# Decide which collections to use
process.IIHEAnalysis.photonCollection       = cms.InputTag('photons'        )
process.IIHEAnalysis.electronCollection     = cms.InputTag('gedGsfElectrons')
#process.IIHEAnalysis.electronCollection     = cms.InputTag('gsfElectrons')
process.IIHEAnalysis.muonCollection         = cms.InputTag('muons'          )
process.IIHEAnalysis.superClusterCollection = cms.InputTag('correctedHybridSuperClusters')

process.IIHEAnalysis.TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')

# Used for the crystal reconstruction in the HEEP module
process.IIHEAnalysis.reducedBarrelRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
process.IIHEAnalysis.reducedEndcapRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')

process.IIHEAnalysis.storeGlobalTrackMuons = cms.untracked.bool(True)
process.IIHEAnalysis.storeStandAloneMuons  = cms.untracked.bool(True)
process.IIHEAnalysis.storeInnerTrackMuons  = cms.untracked.bool(True)

# Triggers:
# Declare whatever lists you like
# Triggers can appear more than once- the analyser is clever enough to only add them once
# You can include wildcard characters to include groups of triggers
HEEPTriggers = []
HEEPTriggers.append('HLT_Ele27_WP80_v8')
HEEPTriggers.append('HLT_DoubleEle33_CaloIdL_v11')
HEEPTriggers.append('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v3')
HEEPTriggers.append('HLT_DoubleEle33_CaloIdT_v7')

# Triggers by topology
singleElectronTriggers = ['singleElectron']
doubleElectronTriggers = ['doubleElectron']
singleMuonTriggers     = ['singleMuon'    ]
doubleMuonTriggers     = ['doubleMuon'    ]
singleElectronSingleMuonTriggers   = ['singleElectronSingleMuon']
singleElectronDoubleMuonTriggers   = ['singleElectronDoubleMuon']
doubleElectronSingleMuonTriggers   = ['doubleElectronSingleMuon']

# Add things together
triggers = singleElectronTriggers + doubleElectronTriggers
#triggers = HEEPTriggers + doubleElectronTriggers

# Now pass the comma separated list of triggers
csvTriggers = ','.join(triggers)
process.IIHEAnalysis.triggers = cms.untracked.string(csvTriggers)

process.IIHEAnalysis.includeEventModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeSuperClusterModule = cms.untracked.bool(True)
process.IIHEAnalysis.includePhotonModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeElectronModule     = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeHEEPModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeMCTruthModule      = cms.untracked.bool(True)
process.IIHEAnalysis.includeTriggerModule      = cms.untracked.bool(False)

process.p1 = cms.Path(process.IIHEAnalysis)

