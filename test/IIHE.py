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
    #'file:/user/aidan/public/Spring14dr__ZPrimePSIToEEMuMu_M-3000_13TeV_pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__067C0233-6ED1-E311-8C07-0025902008EC.root'
    'file:/user/gfasanel/public/0CCBF0FA-2289-E411-A5D4-003048F0E55A.root'
])

# Global tags:
# PHYS14 25ns: PHYS14_25_V1

globalTag = 'PHYS14_25_V1'
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
    input = cms.untracked.int32(100)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('outfile.root')
)

process.load("UserCode.IIHETree.IIHETree_cfi")
process.IIHEAnalysis.globalTag = cms.string(globalTag)
process.IIHEAnalysis.MCTruth_ptThreshold = cms.untracked.double(10.0)
process.IIHEAnalysis.MCTruth_mThreshold  = cms.untracked.double(20.0)
process.IIHEAnalysis.photonCollection    = cms.InputTag('photons'        )
process.IIHEAnalysis.electronCollection  = cms.InputTag('gedGsfElectrons')
process.IIHEAnalysis.muonCollection      = cms.InputTag('muons'          )
process.IIHEAnalysis.reducedBarrelRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
process.IIHEAnalysis.reducedEndcapRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')
process.IIHEAnalysis.muon_triggerDeltaRThreshold = cms.untracked.double(0.5)
process.IIHEAnalysis.HEEP_triggerDeltaRThreshold = cms.untracked.double(0.5)

process.p1 = cms.Path(process.IIHEAnalysis)

