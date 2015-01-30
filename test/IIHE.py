import FWCore.ParameterSet.Config as cms

process = cms.Process("gsfcheckertree")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("RecoTracker.Configuration.RecoTracker_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("CommonTools.ParticleFlow.pfParticleSelection_cff")
process.load("CommonTools.ParticleFlow.Isolation.pfElectronIsolation_cff")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.INFO.limit = 100000
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
    #'file:/user/aidan/public/Spring14dr__ZPrimePSIToEEMuMu_M-3000_13TeV_pythia8__AODSIM__PU20bx25_POSTLS170_V5-v1__067C0233-6ED1-E311-8C07-0025902008EC.root'
    'file:dcap://maite.iihe.ac.be/pnfs/iihe/cms/ph/sc4/store/mc/Fall13dr/ZprimeToMuMu_M-5000_Tune4C_13TeV-pythia8/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/543A5CE4-AF78-E311-9689-0026189438A9.root'
])

# PFMET Type 1 (JEC) correction
process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")

# select global tag and pfJetMET correction automatically from datasetpath when using multicrab
noDataset = True
import sys
for arg in sys.argv:
    if arg.startswith("-CMSSW.datasetpath"):
        dataset = arg[19:]
        noDataset =  False
        break
if noDataset:
    dataset = 'none'

globalTag = 'GR_R_70_V1::All'
process.GlobalTag.globaltag = 'GR_R_70_V1::All'
print "Global Tag is ", process.GlobalTag.globaltag

process.out = cms.OutputModule("PoolOutputModule",
    ##process.FEVTSIMEventContent,
    fileName = cms.untracked.string('outfile.root')
)

process.options = cms.untracked.PSet(
    #fileMode = cms.untracked.string('NOMERGE')
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    #wantSummary = cms.untracked.bool(True)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
    #input = cms.untracked.int32(-1)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('outfile.root')
)

process.hltPhysicsDeclared = cms.EDFilter('HLTPhysicsDeclared',
                                  invert = cms.bool(False),
                                  L1GtReadoutRecordTag = cms.InputTag('gtDigis')
                                  )

## # Primary vertex filter and no scraping events
## # https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCollisionsDataAnalysis
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32 (4),
                                           maxAbsZ = cms.double (24),
                                           maxd0 = cms.double (2)
                                           )
process.primaryVertexPath = cms.Path(process.primaryVertexFilter)

process.noscraping = cms.EDFilter("FilterOutScraping",
                                applyfilter = cms.untracked.bool(True),
                                debugOn = cms.untracked.bool(False),
                                numtrack = cms.untracked.uint32(10),
                                thresh = cms.untracked.double(0.25)
                                )

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.CaloTowerConstituentsMapBuilder = cms.ESProducer("CaloTowerConstituentsMapBuilder",
  MapFile = cms.untracked.string('Geometry/CaloTopology/data/CaloTowerEEGeometric.map.gz')
)

## The next three lines are for rho computation (energy density, highly correlated to PU), see here :
## https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRecipesFor2011#FastJet_based_pile_up_isolation
process.load("RecoJets.JetProducers.kt4PFJets_cfi")
process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

from CommonTools.ParticleFlow.Tools.pfIsolation import setupPFElectronIso, setupPFMuonIso, setupPFPhotonIso
process.eleIsoSequence = cms.Sequence(setupPFElectronIso(process, 'gedGsfElectrons'))
#process.muIsoSequence  = cms.Sequence(    setupPFMuonIso(process, 'muons'          ))
#process.phoIsoSequence = cms.Sequence(  setupPFPhotonIso(process, 'photons'        ))


process.load("UserCode.IIHETree.IIHETree_cfi")
process.IIHEAnalysis.globalTag = cms.string(globalTag)
process.IIHEAnalysis.MCTruth_ptThreshold = cms.untracked.double(10.0)
process.IIHEAnalysis.MCTruth_mThreshold  = cms.untracked.double(20.0)
process.IIHEAnalysis.photonCollection    = cms.InputTag('photons'        )
process.IIHEAnalysis.electronCollection  = cms.InputTag('gedGsfElectrons')
process.IIHEAnalysis.muonCollection      = cms.InputTag('muons'          )

process.otherStuff = cms.Sequence( process.kt6PFJets )

process.load("RecoMET.METFilters.ecalLaserCorrFilter_cfi")
process.load('RecoMET.METFilters.eeBadScFilter_cfi')
process.MessageLogger.suppressError = cms.untracked.vstring ('ecalLaserCorrFilter')

process.p1 = cms.Path(process.otherStuff * process.hltPhysicsDeclared * process.eeBadScFilter * process.ecalLaserCorrFilter * process.noscraping * process.primaryVertexFilter * process.pfParticleSelectionSequence * process.eleIsoSequence * process.producePFMETCorrections * process.IIHEAnalysis)


