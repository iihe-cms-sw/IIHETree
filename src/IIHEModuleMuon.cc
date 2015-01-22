#include "UserCode/IIHETree/interface/IIHEModuleMuon.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleMuon::IIHEModuleMuon(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
  triggerDeltaRThreshold_ = iConfig.getUntrackedParameter<double>("muon_triggerDeltaRThreshold", 1.0) ;
}
IIHEModuleMuon::~IIHEModuleMuon(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleMuon::beginJob(){
  addBranch("muon_n", kUInt) ;
  addBranch("muon_charge", kVectorInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("muon_qoverp") ;
  addBranch("muon_pt") ;
  addBranch("muon_eta") ;
  addBranch("muon_phi") ;
  addBranch("muon_p") ;
  addBranch("muon_px") ;
  addBranch("muon_py") ;
  addBranch("muon_pz") ;
  addBranch("muon_theta") ;
  addBranch("muon_lambda") ;
  addBranch("muon_dxy") ;
  addBranch("muon_d0" ) ;
  addBranch("muon_dsz") ;
  addBranch("muon_dz" ) ;
  addBranch("muon_dxy_beamSpot" ) ;
  addBranch("muon_dxy_firstPVtx" ) ;
  addBranch("muon_dz_beamSpot" ) ;
  addBranch("muon_dz_firstPVtx" ) ;
  addBranch("muon_vx" ) ;
  addBranch("muon_vy" ) ;
  addBranch("muon_vz" ) ;
  
  addBranch("muon_outerPt") ;
  addBranch("muon_outerEta") ;
  addBranch("muon_outerPhi") ;
  addBranch("muon_outerTheta") ;
  setBranchType(kVectorInt) ;
  addBranch("muon_numberOfValidPixelHits") ;
  addBranch("muon_numberOfValidTrackerHits") ;
  addBranch("muon_numberOfValidMuonHits") ;
  addBranch("muon_numberOfValidHits") ;
  addBranch("muon_trackerLayersWithMeasurement") ;
  addBranch("muon_numberOfLostHits") ;
  addBranch("muon_numberOfMatchedStations") ;
  addBranch("muon_validFraction", kVectorDouble) ;
  
  setBranchType(kVectorBool) ;
  addBranch("muon_isGlobalMuon") ;
  addBranch("muon_isTrackerMuon") ;
  addBranch("muon_isPFMuon") ;
  addBranch("muon_isPFIsolationValid") ;
  
  setBranchType(kVectorFloat) ;
  addBranch("muon_chi2") ;
  addBranch("muon_ndof") ;
  addBranch("muon_normalizedChi2") ;
  
  addBranch("muon_qoverpError") ;
  addBranch("muon_ptError") ;
  addBranch("muon_thetaError") ;
  addBranch("muon_lambdaError") ;
  addBranch("muon_phiError") ;
  addBranch("muon_dxyError") ;
  addBranch("muon_d0Error") ;
  addBranch("muon_dszError") ;
  addBranch("muon_dzError") ;
  addBranch("muon_etaError") ;
  
  // Isolation variables
  setBranchType(kVectorFloat) ;
  addBranch("muon_isolationR03_sumPt") ;
  addBranch("muon_isolationR03_trackerVetoPt") ;
  addBranch("muon_isolationR03_emEt") ;
  addBranch("muon_isolationR03_emVetoEt") ;
  addBranch("muon_isolationR03_hadEt") ;
  addBranch("muon_isolationR03_hadVetoEt") ;
  
  addBranch("muon_isolationR05_sumPt") ;
  addBranch("muon_isolationR05_trackerVetoPt") ;
  addBranch("muon_isolationR05_emEt") ;
  addBranch("muon_isolationR05_emVetoEt") ;
  addBranch("muon_isolationR05_hadEt") ;
  addBranch("muon_isolationR05_hadVetoEt") ;
  
  addBranch("muon_pfIsolationR03_sumChargedHadronPt"             ) ;
  addBranch("muon_pfIsolationR03_sumChargedParticlePt"           ) ;
  addBranch("muon_pfIsolationR03_sumPhotonEt"                    ) ;
  addBranch("muon_pfIsolationR03_sumNeutralHadronEtHighThreshold") ;
  addBranch("muon_pfIsolationR03_sumPhotonEtHighThreshold"       ) ;
  addBranch("muon_pfIsolationR03_sumPUPt"                        ) ;
  
  addBranch("muon_pfIsolationR04_sumChargedHadronPt"             ) ;
  addBranch("muon_pfIsolationR04_sumChargedParticlePt"           ) ;
  addBranch("muon_pfIsolationR04_sumPhotonEt"                    ) ;
  addBranch("muon_pfIsolationR04_sumNeutralHadronEtHighThreshold") ;
  addBranch("muon_pfIsolationR04_sumPhotonEtHighThreshold"       ) ;
  addBranch("muon_pfIsolationR04_sumPUPt"                        ) ;
  
  addBranch("muon_pfMeanDRIsoProfileR03_sumChargedHadronPt"             ) ;
  addBranch("muon_pfMeanDRIsoProfileR03_sumChargedParticlePt"           ) ;
  addBranch("muon_pfMeanDRIsoProfileR03_sumPhotonEt"                    ) ;
  addBranch("muon_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold") ;
  addBranch("muon_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold"       ) ;
  addBranch("muon_pfMeanDRIsoProfileR03_sumPUPt"                        ) ;
  
  addBranch("muon_pfMeanDRIsoProfileR04_sumChargedHadronPt"             ) ;
  addBranch("muon_pfMeanDRIsoProfileR04_sumChargedParticlePt"           ) ;
  addBranch("muon_pfMeanDRIsoProfileR04_sumPhotonEt"                    ) ;
  addBranch("muon_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold") ;
  addBranch("muon_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold"       ) ;
  addBranch("muon_pfMeanDRIsoProfileR04_sumPUPt"                        ) ;
  
  addBranch("muon_innerPosition_x") ;
  addBranch("muon_innerPosition_y") ;
  addBranch("muon_innerPosition_z") ;
  
  // TeV optimized values
  addBranch("muon_tevOptimized_charge", kVectorInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("muon_tevOptimized_pt") ;
  addBranch("muon_tevOptimized_eta") ;
  addBranch("muon_tevOptimized_phi") ;
  addBranch("muon_tevOptimized_theta") ;
  addBranch("muon_tevOptimized_px") ;
  addBranch("muon_tevOptimized_py") ;
  addBranch("muon_tevOptimized_pz") ;
  addBranch("muon_tevOptimized_d0") ;
  addBranch("muon_tevOptimized_dz") ;
  addBranch("muon_tevOptimized_dz_beamSpot") ;
  addBranch("muon_tevOptimized_dz_firstPVtx") ;
  addBranch("muon_tevOptimized_dxy") ;
  addBranch("muon_tevOptimized_dxy_beamSpot") ;
  addBranch("muon_tevOptimized_dxy_firstPVtx") ;
  
  addBranch("muon_tevOptimized_ptError") ;
  addBranch("muon_tevOptimized_etaError") ;
  addBranch("muon_tevOptimized_phiError") ;
  addBranch("muon_tevOptimized_thetaError") ;
  addBranch("muon_tevOptimized_d0Error") ;
  addBranch("muon_tevOptimized_dzError") ;
  addBranch("muon_tevOptimized_dxyError") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                      Triggers                                      //
  ////////////////////////////////////////////////////////////////////////////////////////
  addTriggerHLTMuon("hltL1sMu16Eta2p1"                          , triggerDeltaRThreshold_) ;
  addTriggerHLTMuon("hltL1sL1Mu3p5EG12"                         , triggerDeltaRThreshold_) ;
  addTriggerHLTMuon("hltL1Mu3p5EG12L3Filtered22"                , triggerDeltaRThreshold_) ;
  addTriggerHLTMuon("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q", triggerDeltaRThreshold_) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Get the beamspot from the Event:
  // The beamspot is passed to the IIHEAnalysis class, so we call it from parent_
  // Don't forget to declare IIHEModuleVertex as a friend of IIHEAnalysis!
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(parent_->beamSpotLabel_, theBeamSpot);

  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint* firstPrimaryVertex = parent_->getFirstPrimaryVertex() ;

  // Trigger information
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag,trigEvent);
  
  // Muon collections
  reco::MuonCollection muons = parent_->getMuonCollection() ;
  
  store("muon_n", (unsigned int)(muons.size())) ;
  for(reco::MuonCollection::const_iterator muIt = muons.begin(); muIt != muons.end(); ++muIt) {
    if(muIt->isGlobalMuon() || true){
      // get TeV optimized track
      reco::Muon::MuonTrackTypePair tevOptimizedTrk = muon::tevOptimized(*muIt, 200, 17., 40., 0.25) ;
      
      store("muon_charge"       , muIt->globalTrack()->charge()                ) ;
      store("muon_qoverp"       , muIt->globalTrack()->qoverp()                ) ;
      store("muon_pt"           , muIt->globalTrack()->pt()                    ) ;
      store("muon_eta"          , muIt->globalTrack()->eta()                   ) ;
      store("muon_phi"          , muIt->globalTrack()->phi()                   ) ;
      store("muon_p"            , muIt->globalTrack()->p()                     ) ;
      store("muon_px"           , muIt->globalTrack()->px()                    ) ;
      store("muon_py"           , muIt->globalTrack()->py()                    ) ;
      store("muon_pz"           , muIt->globalTrack()->pz()                    ) ;
      store("muon_theta"        , muIt->globalTrack()->theta()                 ) ;
      store("muon_lambda"       , muIt->globalTrack()->lambda()                ) ;
      store("muon_d0"           , muIt->globalTrack()->d0()                    ) ;
      store("muon_dz"           , muIt->globalTrack()->dz()                    ) ;
      store("muon_dz_beamSpot"  , muIt->globalTrack()->dz(beamspot)            ) ;
      store("muon_dz_firstPVtx" , muIt->globalTrack()->dz(*firstPrimaryVertex) ) ;
      store("muon_dxy"          , muIt->globalTrack()->dxy()                   ) ;
      store("muon_dxy_beamSpot" , muIt->globalTrack()->dxy(beamspot)           ) ;
      store("muon_dxy_firstPVtx", muIt->globalTrack()->dxy(*firstPrimaryVertex)) ;
      store("muon_dxy"          , muIt->globalTrack()->dxy(beamspot)           ) ;
      store("muon_dsz"          , muIt->globalTrack()->dsz(beamspot)           ) ;
      store("muon_dz"           , muIt->globalTrack()->dz(beamspot)            ) ;
      store("muon_vx"           , muIt->globalTrack()->vx()                    ) ;
      store("muon_vy"           , muIt->globalTrack()->vy()                    ) ;
      store("muon_vz"           , muIt->globalTrack()->vz()                    ) ;
      store("muon_outerPt"      , muIt->globalTrack()->outerPt()               ) ;
      store("muon_outerEta"     , muIt->globalTrack()->outerEta()              ) ;
      store("muon_outerPhi"     , muIt->globalTrack()->outerPhi()              ) ;
      store("muon_outerTheta"   , muIt->globalTrack()->outerTheta()            ) ;
      
      float etaError = muIt->globalTrack()->thetaError()/sin(muIt->globalTrack()->theta()) ;
      store("muon_qoverpError", muIt->globalTrack()->qoverpError()) ;
      store("muon_ptError"    , muIt->globalTrack()->ptError()    ) ;
      store("muon_thetaError" , muIt->globalTrack()->thetaError() ) ;
      store("muon_lambdaError", muIt->globalTrack()->lambdaError()) ;
      store("muon_phiError"   , muIt->globalTrack()->phiError()   ) ;
      store("muon_dxyError"   , muIt->globalTrack()->dxyError()   ) ;
      store("muon_d0Error"    , muIt->globalTrack()->d0Error()    ) ;
      store("muon_dszError"   , muIt->globalTrack()->dszError()   ) ;
      store("muon_dzError"    , muIt->globalTrack()->dzError()    ) ;
      store("muon_etaError"   , etaError                          ) ;
      
      store("muon_numberOfValidPixelHits"      , muIt->innerTrack()->hitPattern().numberOfValidPixelHits()   ) ;
      store("muon_numberOfValidTrackerHits"    , muIt->globalTrack()->hitPattern().numberOfValidTrackerHits()) ;
      store("muon_numberOfValidMuonHits"       , muIt->globalTrack()->hitPattern().numberOfValidMuonHits()   ) ;
      store("muon_numberOfValidHits"           , muIt->globalTrack()->numberOfValidHits()                    ) ;
      store("muon_trackerLayersWithMeasurement", muIt->track()->hitPattern().trackerLayersWithMeasurement()  ) ;
      store("muon_numberOfLostHits"            , muIt->globalTrack()->numberOfLostHits()                     ) ;
      store("muon_numberOfMatchedStations"     , muIt->numberOfMatchedStations()                             ) ;
      store("muon_validFraction"               , muIt->globalTrack()->validFraction()                        ) ;
             
      store("muon_isGlobalMuon"      , muIt->isGlobalMuon()                    ) ;        
      store("muon_isTrackerMuon"     , muIt->isTrackerMuon()                   ) ;        
      store("muon_isPFMuon"          , muIt->isPFMuon()                        ) ;        
      store("muon_isPFIsolationValid", muIt->isPFIsolationValid()              ) ;        
      store("muon_chi2"              , muIt->globalTrack()->chi2()             ) ;
      store("muon_ndof"              , muIt->globalTrack()->ndof()             ) ;
      store("muon_normalizedChi2"    , muIt->globalTrack()->normalizedChi2()   ) ;
      store("muon_innerPosition_x"   , muIt->globalTrack()->innerPosition().X()) ;
      store("muon_innerPosition_y"   , muIt->globalTrack()->innerPosition().Y()) ;
      store("muon_innerPosition_z"   , muIt->globalTrack()->innerPosition().Z()) ;
      
      store("muon_isolationR03_sumPt"        , muIt->isolationR03().sumPt              ) ;
      store("muon_isolationR03_trackerVetoPt", muIt->isolationR03().trackerVetoPt      ) ;
      store("muon_isolationR03_emEt"         , muIt->isolationR03().emEt               ) ;
      store("muon_isolationR03_emVetoEt"     , muIt->isolationR03().emVetoEt           ) ;
      store("muon_isolationR03_hadEt"        , muIt->isolationR03().hadEt              ) ;
      store("muon_isolationR03_hadVetoEt"    , muIt->isolationR03().hadVetoEt          ) ;
      
      store("muon_isolationR05_sumPt"        , muIt->isolationR05().sumPt              ) ;
      store("muon_isolationR05_trackerVetoPt", muIt->isolationR05().trackerVetoPt      ) ;
      store("muon_isolationR05_emEt"         , muIt->isolationR05().emEt               ) ;
      store("muon_isolationR05_emVetoEt"     , muIt->isolationR05().emVetoEt           ) ;
      store("muon_isolationR05_hadEt"        , muIt->isolationR05().hadEt              ) ;
      store("muon_isolationR05_hadVetoEt"    , muIt->isolationR05().hadVetoEt          ) ;
      
      store("muon_pfIsolationR03_sumChargedHadronPt"             , muIt->pfIsolationR03().sumChargedHadronPt             ) ;
      store("muon_pfIsolationR03_sumChargedParticlePt"           , muIt->pfIsolationR03().sumChargedParticlePt           ) ;
      store("muon_pfIsolationR03_sumPhotonEt"                    , muIt->pfIsolationR03().sumPhotonEt                    ) ;
      store("muon_pfIsolationR03_sumNeutralHadronEtHighThreshold", muIt->pfIsolationR03().sumNeutralHadronEtHighThreshold) ;
      store("muon_pfIsolationR03_sumPhotonEtHighThreshold"       , muIt->pfIsolationR03().sumPhotonEtHighThreshold       ) ;
      store("muon_pfIsolationR03_sumPUPt"                        , muIt->pfIsolationR03().sumPUPt                        ) ;
      
      store("muon_pfIsolationR04_sumChargedHadronPt"             , muIt->pfIsolationR03().sumChargedHadronPt             ) ;
      store("muon_pfIsolationR04_sumChargedParticlePt"           , muIt->pfIsolationR03().sumChargedParticlePt           ) ;
      store("muon_pfIsolationR04_sumPhotonEt"                    , muIt->pfIsolationR03().sumPhotonEt                    ) ;
      store("muon_pfIsolationR04_sumNeutralHadronEtHighThreshold", muIt->pfIsolationR03().sumNeutralHadronEtHighThreshold) ;
      store("muon_pfIsolationR04_sumPhotonEtHighThreshold"       , muIt->pfIsolationR03().sumPhotonEtHighThreshold       ) ;
      store("muon_pfIsolationR04_sumPUPt"                        , muIt->pfIsolationR03().sumPUPt                        ) ;
      
      store("muon_pfMeanDRIsoProfileR03_sumChargedHadronPt"             , muIt->pfMeanDRIsoProfileR03().sumChargedHadronPt             ) ;
      store("muon_pfMeanDRIsoProfileR03_sumChargedParticlePt"           , muIt->pfMeanDRIsoProfileR03().sumChargedParticlePt           ) ;
      store("muon_pfMeanDRIsoProfileR03_sumPhotonEt"                    , muIt->pfMeanDRIsoProfileR03().sumPhotonEt                    ) ;
      store("muon_pfMeanDRIsoProfileR03_sumNeutralHadronEtHighThreshold", muIt->pfMeanDRIsoProfileR03().sumNeutralHadronEtHighThreshold) ;
      store("muon_pfMeanDRIsoProfileR03_sumPhotonEtHighThreshold"       , muIt->pfMeanDRIsoProfileR03().sumPhotonEtHighThreshold       ) ;
      store("muon_pfMeanDRIsoProfileR03_sumPUPt"                        , muIt->pfMeanDRIsoProfileR03().sumPUPt                        ) ;
      
      store("muon_pfMeanDRIsoProfileR04_sumChargedHadronPt"             , muIt->pfMeanDRIsoProfileR04().sumChargedHadronPt             ) ;
      store("muon_pfMeanDRIsoProfileR04_sumChargedParticlePt"           , muIt->pfMeanDRIsoProfileR04().sumChargedParticlePt           ) ;
      store("muon_pfMeanDRIsoProfileR04_sumPhotonEt"                    , muIt->pfMeanDRIsoProfileR04().sumPhotonEt                    ) ;
      store("muon_pfMeanDRIsoProfileR04_sumNeutralHadronEtHighThreshold", muIt->pfMeanDRIsoProfileR04().sumNeutralHadronEtHighThreshold) ;
      store("muon_pfMeanDRIsoProfileR04_sumPhotonEtHighThreshold"       , muIt->pfMeanDRIsoProfileR04().sumPhotonEtHighThreshold       ) ;
      store("muon_pfMeanDRIsoProfileR04_sumPUPt"                        , muIt->pfMeanDRIsoProfileR04().sumPUPt                        ) ;
      
      store("muon_tevOptimized_charge"       , tevOptimizedTrk.first->charge()                ) ;
      store("muon_tevOptimized_pt"           , tevOptimizedTrk.first->pt()                    ) ;
      store("muon_tevOptimized_eta"          , tevOptimizedTrk.first->eta()                   ) ;
      store("muon_tevOptimized_phi"          , tevOptimizedTrk.first->phi()                   ) ;
      store("muon_tevOptimized_theta"        , tevOptimizedTrk.first->theta()                 ) ;
      store("muon_tevOptimized_px"           , tevOptimizedTrk.first->px()                    ) ;
      store("muon_tevOptimized_py"           , tevOptimizedTrk.first->py()                    ) ;
      store("muon_tevOptimized_pz"           , tevOptimizedTrk.first->pz()                    ) ;
      store("muon_tevOptimized_d0"           , tevOptimizedTrk.first->d0()                    ) ;
      store("muon_tevOptimized_dz"           , tevOptimizedTrk.first->dz()                    ) ;
      store("muon_tevOptimized_dz_beamSpot"  , tevOptimizedTrk.first->dz(beamspot)            ) ;
      store("muon_tevOptimized_dz_firstPVtx" , tevOptimizedTrk.first->dz(*firstPrimaryVertex) ) ;
      store("muon_tevOptimized_dxy"          , tevOptimizedTrk.first->dxy()                   ) ;
      store("muon_tevOptimized_dxy_beamSpot" , tevOptimizedTrk.first->dxy(beamspot)           ) ;
      store("muon_tevOptimized_dxy_firstPVtx", tevOptimizedTrk.first->dxy(*firstPrimaryVertex)) ;
      
      store("muon_tevOptimized_ptError"   , tevOptimizedTrk.first->ptError()   ) ;
      store("muon_tevOptimized_etaError"  , tevOptimizedTrk.first->etaError()  ) ;
      store("muon_tevOptimized_phiError"  , tevOptimizedTrk.first->phiError()  ) ;
      store("muon_tevOptimized_thetaError", tevOptimizedTrk.first->thetaError()) ;
      store("muon_tevOptimized_d0Error"   , tevOptimizedTrk.first->d0Error()   ) ;
      store("muon_tevOptimized_dzError"   , tevOptimizedTrk.first->dzError()   ) ;
      store("muon_tevOptimized_dxyError"  , tevOptimizedTrk.first->dxyError()  ) ;
    }
  }
}

void IIHEModuleMuon::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMuon::beginEvent(){}
void IIHEModuleMuon::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMuon::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMuon);
