##########################################################################################
#                                      Global tags                                       #
##########################################################################################
# Global tags:
# PHYS14 25ns: PHYS14_25_V1

#globalTag = 'MCRUN2_74_V9::All'  # 25ns asymptotic
#globalTag = 'MCRUN2_74_V9A::All' # 50ns asymptotic
#globalTag = '74X_mcRun2_asymptotic_realisticBS_v1' # 25ns realistic beamspot
#globalTag = '741_p1_mcRun2_Realistic_50ns_v0'      # 50ns realistic beamspot

#globalTag = 'PHYS14_25_V1::All'
#globalTag = 'GR_R_70_V1::All'

# Data
#globalTag = 'GR_H_V58C::All' # 2015, AB
#globalTag = '74X_dataRun2_HLT_v0' # 2015, C
#
#globalTag='MCRUN2_74_V9A::All'
#globalTag = '74X_dataRun2_Prompt_v1' # Run2015B, Run2015C
globalTag = '74X_dataRun2_Prompt_v2' # Run2015D


##########################################################################################
#                                  Start the sequences                                   #
##########################################################################################
import FWCore.ParameterSet.Config as cms

process = cms.Process("IIHEAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.EventContent.EventContent_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100

##########################################################################################
#                                         Files                                          #
##########################################################################################
readFiles = cms.untracked.vstring()
secFiles  = cms.untracked.vstring()
process.source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend([
'/store/data/Run2015D/DoubleEG/AOD/PromptReco-v3/000/256/630/00000/C4A65920-395F-E511-B511-02163E01264C.root'
])

filename_out = 'outfile.root'
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string(filename_out) )
process.TFileService = cms.Service("TFileService", fileName = cms.string(          filename_out) )

##########################################################################################
#                                     Main options                                       #
##########################################################################################
process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag

process.options = cms.untracked.PSet( SkipEvent = cms.untracked.vstring('ProductNotFound') )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Taken from Sherif's code
from RecoJets.JetProducers.kt4PFJets_cfi import *
process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

##########################################################################################
#                                   IIHETree options                                     #
##########################################################################################
process.load("UserCode.IIHETree.IIHETree_cfi")

# Only save some triggers.
process.IIHEAnalysis.TriggerResults = cms.InputTag('TriggerResults', '', 'HLT')
triggers = 'singleElectron;doubleElectron;singleMuon;singleElectronSingleMuon;singleElectronDoubleMuon;doubleElectronSingleMuon'
process.IIHEAnalysis.triggers = cms.untracked.string(triggers)

#process.IIHEAnalysis.triggers = cms.untracked.string('doubleElectron')

process.IIHEAnalysis.globalTag = cms.string(globalTag)

process.IIHEAnalysis.eventRho = cms.InputTag('kt6PFJetsForIsolation:rho')

# Collections.
process.IIHEAnalysis.photonCollection    = cms.InputTag('photons'        )
process.IIHEAnalysis.electronCollection  = cms.InputTag('gedGsfElectrons')
process.IIHEAnalysis.muonCollection      = cms.InputTag('muons'          )
process.IIHEAnalysis.superClusterCollection = cms.InputTag('correctedHybridSuperClusters')
process.IIHEAnalysis.reducedBarrelRecHitCollection = cms.InputTag('reducedEcalRecHitsEB')
process.IIHEAnalysis.reducedEndcapRecHitCollection = cms.InputTag('reducedEcalRecHitsEE')

miniAOD = False
if miniAOD:
    process.IIHEAnalysis.primaryVertex = cms.InputTag('offlineSlimmedPrimaryVertices')
    process.IIHEAnalysis.photonCollection    = cms.InputTag('slimmedPhotons'  )
    process.IIHEAnalysis.electronCollection  = cms.InputTag('slimmedElectrons')
    process.IIHEAnalysis.muonCollection      = cms.InputTag('slimmedMuons'    )

# Trigger matching stuff.  0.5 should be sufficient.
process.IIHEAnalysis.muon_triggerDeltaRThreshold = cms.untracked.double(0.5)
process.IIHEAnalysis.HEEP_triggerDeltaRThreshold = cms.untracked.double(0.5)

# In the absence of high ET electrons, only save events with really high Z candidates.
process.IIHEAnalysis.ZBosonZMassAcceptLower    = cms.untracked.double(850)
# Don't bother with J/psi or Upsilon, they will only weigh us down!
process.IIHEAnalysis.ZBosonJPsiAcceptMassLower = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonJPsiAcceptMassUpper = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonUpsAcceptMassLower  = cms.untracked.double(1e6)
process.IIHEAnalysis.ZBosonUpsAcceptMassUpper  = cms.untracked.double(1e6)

# But make sure we save Z bosons from 50 GeV and up.
process.IIHEAnalysis.ZBosonZMassLowerCuttoff   = cms.untracked.double( 50)

process.IIHEAnalysis.ZBosonDeltaRCut           = cms.untracked.double(1e-3)

# Only save Z->ee, Z->em.
pt_threshold = 5
process.IIHEAnalysis.ZBosonEtThreshold = cms.untracked.double(pt_threshold)
process.IIHEAnalysis.ZBosonSaveZee  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZmm  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZem  = cms.untracked.bool(True )
process.IIHEAnalysis.ZBosonSaveZeeg = cms.untracked.bool(False)
process.IIHEAnalysis.ZBosonSaveZmmg = cms.untracked.bool(False)

process.IIHEAnalysis.electrons_ETThreshold = cms.untracked.double(pt_threshold)
process.IIHEAnalysis.muon_pTThreshold      = cms.untracked.double(pt_threshold)

process.IIHEAnalysis.LeptonsAccept_pTThreshold = cms.untracked.double(10)
# Require at least two leptons...
process.IIHEAnalysis.LeptonsAccept_nLeptons    = cms.untracked.double(2)
# ...at least one of which is an electron.
process.IIHEAnalysis.LeptonsAccept_nElectrons  = cms.untracked.double(1)

process.IIHEAnalysis.includeTriggerModule         = cms.untracked.bool(True )
# Turn off these guys, since we're working with data.
process.IIHEAnalysis.includeMCTruthModule         = cms.untracked.bool(False)
process.IIHEAnalysis.includeAutoAcceptEventModule = cms.untracked.bool(False)

process.IIHEAnalysis.debug = cms.bool(False)

##########################################################################################
#                            Woohoo!  We're ready to start!                              #
##########################################################################################
process.p1 = cms.Path(process.kt6PFJetsForIsolation+process.IIHEAnalysis)

