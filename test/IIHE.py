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
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/0AD36E59-BD6B-E411-BABC-00266CFFA7A8.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/16D0A2E6-BB6B-E411-B3AF-00266CF32920.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/4445BF5A-A66B-E411-B462-002590AC4E28.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/9AB68F75-396C-E411-A28D-002590DB9262.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/00000/A8367075-396C-E411-8FE6-002590DB9262.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/46FE2F45-BD6B-E411-83B5-008CFA104E64.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/6628077B-016C-E411-A2B5-00266CFFA754.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/7C413593-BE6B-E411-B8E6-0025904B1420.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/82BD2B28-B16B-E411-9EB8-002481E94C56.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/9085294F-016C-E411-9A4C-0025907DCA9C.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/9AA8B9D8-D36B-E411-AABA-D8D385FF6C5E.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/BAE5B92F-B66B-E411-B7B6-002481E101DC.root',
    '/store/mc/Phys14DR/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/AODSIM/PU20bx25_tsg_PHYS14_25_V1-v1/10000/FA016B56-AA6B-E411-8D1A-00266CF25E44.root'
    #'file:/user/aidan/public/Spring14dr__ZPrimePSIToEEMuMu_M-3000_13TeV_pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__067C0233-6ED1-E311-8C07-0025902008EC.root'
    #'file:/user/gfasanel/public/0CCBF0FA-2289-E411-A5D4-003048F0E55A.root'
    #'file:Run2012A.root'
    #'file:SingleElectron.root'
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
process.IIHEAnalysis.includeTriggerModule      = cms.untracked.bool(True)

process.p1 = cms.Path(process.IIHEAnalysis)

