import FWCore.ParameterSet.Config as cms

IIHEAnalysis = cms.EDAnalyzer("IIHEAnalysis",
    debug    = cms.bool(True),
    beamSpot = cms.InputTag("offlineBeamSpot")
)

