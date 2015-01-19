import FWCore.ParameterSet.Config as cms
import getpass, os

pwd = os.getcwd()
os.chdir("UserCode/IIHETree/")
username = getpass.getuser()
os.system("git log -n 1 | head -n 1 | awk '{print $2}' > /tmp/%s_git.hash"%username)
f = open('/tmp/git.hash')
git_hash = f.read().rstrip('\n')
print 'Using git hash: ' , git_hash
os.chdir(pwd)

IIHEAnalysis = cms.EDAnalyzer("IIHEAnalysis",
    debug         = cms.bool(True),
    beamSpot      = cms.InputTag("offlineBeamSpot"),
    primaryVertex = cms.InputTag('offlinePrimaryVertices'),
    git_hash = cms.string(git_hash),
    EcalHcal1EffAreaBarrel  = cms.untracked.double(0.28),
    EcalHcal1EffAreaEndcaps = cms.untracked.double(0.28)
)

