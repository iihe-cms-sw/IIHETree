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
    'file:/user/aidan/public/Spring14dr__ZPrimePSIToEEMuMu_M-3000_13TeV_pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__067C0233-6ED1-E311-8C07-0025902008EC.root'
    #'file:/user/gfasanel/public/0CCBF0FA-2289-E411-A5D4-003048F0E55A.root'
])

# Global tags:
# PHYS14 25ns: PHYS14_25_V1

globalTag = 'PHYS14_25_V1::All'
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
    input = cms.untracked.int32(1000)
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
process.IIHEAnalysis.photonCollection   = cms.InputTag('photons'        )
process.IIHEAnalysis.electronCollection = cms.InputTag('gedGsfElectrons')
process.IIHEAnalysis.muonCollection     = cms.InputTag('muons'          )

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
HEEPTriggers.append('HLT_Ele27_WP80_v13')
HEEPTriggers.append('HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v1')

# Triggers by topology
singleElectronTriggers = ['singleElectron']
doubleElectronTriggers = ['doubleElectron']
singleMuonTriggers     = ['singleMuon'    ]
doubleMuonTriggers     = ['doubleMuon'    ]
singleElectronSingleMuonTriggers   = ['singleElectronSingleMuon']
singleElectronDoubleMuonTriggers   = ['singleElectronDoubleMuon']
doubleElectronSingleMuonTriggers   = ['soubleElectronSingleMuon']

# Add things together
triggers = HEEPTriggers + singleElectronTriggers + doubleElectronTriggers

triggers = HEEPTriggers + doubleElectronTriggers

# Now pass the comma separated list of triggers
csvTriggers = ','.join(triggers)
process.IIHEAnalysis.triggers = cms.untracked.string(csvTriggers)

# For trigger matching you need to pass the list of triggers you want to match
# At the moment I'm limiting things to just two triggers
csvElectronTriggerMatching = ','.join(HEEPTriggers)
process.IIHEAnalysis.triggerPhotonMatchings   = cms.untracked.string(csvElectronTriggerMatching)
process.IIHEAnalysis.triggerElectronMatchings = cms.untracked.string(csvElectronTriggerMatching)
process.IIHEAnalysis.triggerMuonMatchings     = cms.untracked.string(csvElectronTriggerMatching)
process.IIHEAnalysis.triggerTauMatchings      = cms.untracked.string(csvElectronTriggerMatching)
process.IIHEAnalysis.triggerJetMatchings      = cms.untracked.string(csvElectronTriggerMatching)

process.IIHEAnalysis.includeEventModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeSuperClusterModule = cms.untracked.bool(True)
process.IIHEAnalysis.includePhotonModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeElectronModule     = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeHEEPModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeMCTruthModule      = cms.untracked.bool(True)
process.IIHEAnalysis.includeTriggerModule      = cms.untracked.bool(True)

process.p1 = cms.Path(process.IIHEAnalysis)

