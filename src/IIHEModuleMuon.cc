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

IIHEModuleMuon::IIHEModuleMuon(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
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
  addBranch("muon_dxy_firstPVtxwithBS" ) ;
  addBranch("muon_dz_beamSpot" ) ;
  addBranch("muon_dz_firstPVtx" ) ;
  addBranch("muon_dz_firstPVtxwithBS" ) ;
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
  addBranch("muon_tevOptimized_dz_firstPVtxwithBS") ;
  addBranch("muon_tevOptimized_dxy") ;
  addBranch("muon_tevOptimized_dxy_beamSpot") ;
  addBranch("muon_tevOptimized_dxy_firstPVtx") ;
  addBranch("muon_tevOptimized_dxy_firstPVtxwithBS") ;
  
  addBranch("muon_tevOptimized_ptError") ;
  addBranch("muon_tevOptimized_etaError") ;
  addBranch("muon_tevOptimized_phiError") ;
  addBranch("muon_tevOptimized_thetaError") ;
  addBranch("muon_tevOptimized_d0Error") ;
  addBranch("muon_tevOptimized_dzError") ;
  addBranch("muon_tevOptimized_dxyError") ;
  
  setBranchType(kVectorBool) ;
  addBranch("muMatch_hltL1sMu16Eta2p1") ;
  addBranch("muMatch_hltL1sL1Mu3p5EG12") ;
  addBranch("muMatch_hltL1Mu3p5EG12L3Filtered22") ;
  addBranch("muMatch_hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q") ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Get the beamspot from the Event:
  // The beamspot is passed to the IIHEAnalysis class, so we call it from parent_
  // Don't forget to declare IIHEModuleVertex as a friend of IIHEAnalysis!
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(parent_->beamSpotLabel_, theBeamSpot);

  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint firstpvertex      (0.0,0.0,0.0) ;
  math::XYZPoint firstpvertexwithBS(0.0,0.0,0.0) ;

  // Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();
  
  Handle<reco::VertexCollection> primaryVertexCollwithBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVertexCollwithBS);
  const reco::VertexCollection* pvcollwithBS = primaryVertexCollwithBS.product();
  
  if(pvcoll->size()>0){
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  if(pvcollwithBS->size()>0){
    reco::VertexCollection::const_iterator firstpv = pvcollwithBS->begin();
    firstpvertexwithBS.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  // Trigger information
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag,trigEvent);
  
  
  // Muon collections
  reco::MuonCollection muons = parent_->getMuonCollection() ;
  
  store("muon_n", (unsigned int)(muons.size())) ;
  for(reco::MuonCollection::const_iterator muIt = muons.begin(); muIt != muons.end(); ++muIt) {
    if (muIt->isGlobalMuon()) {
      // get TeV optimized track
      reco::Muon::MuonTrackTypePair tevOptimizedTrk = muon::tevOptimized(*muIt, 200, 17., 40., 0.25) ;
      
      store("muon_charge"                          , muIt->globalTrack()->charge()                                    ) ;
      store("muon_qoverp"                          , muIt->globalTrack()->qoverp()                                    ) ;
      store("muon_pt"                              , muIt->globalTrack()->pt()                                        ) ;
      store("muon_eta"                             , muIt->globalTrack()->eta()                                       ) ;
      store("muon_phi"                             , muIt->globalTrack()->phi()                                       ) ;
      store("muon_p"                               , muIt->globalTrack()->p()                                         ) ;
      store("muon_px"                              , muIt->globalTrack()->px()                                        ) ;
      store("muon_py"                              , muIt->globalTrack()->py()                                        ) ;
      store("muon_pz"                              , muIt->globalTrack()->pz()                                        ) ;
      store("muon_theta"                           , muIt->globalTrack()->theta()                                     ) ;
      store("muon_lambda"                          , muIt->globalTrack()->lambda()                                    ) ;
      store("muon_d0"                              , muIt->globalTrack()->d0()                                        ) ;
      store("muon_dz"                              , muIt->globalTrack()->dz()                                        ) ;
      store("muon_dz_beamSpot"                     , muIt->globalTrack()->dz(beamspot)                                ) ;
      store("muon_dz_firstPVtx"                    , muIt->globalTrack()->dz(firstpvertex)                            ) ;
      store("muon_dz_firstPVtxwithBS"              , muIt->globalTrack()->dz(firstpvertexwithBS)                      ) ;
      store("muon_dxy"                             , muIt->globalTrack()->dxy()                                       ) ;
      store("muon_dxy_beamSpot"                    , muIt->globalTrack()->dxy(beamspot)                               ) ;
      store("muon_dxy_firstPVtx"                   , muIt->globalTrack()->dxy(firstpvertex)                           ) ;
      store("muon_dxy_firstPVtxwithBS"             , muIt->globalTrack()->dxy(firstpvertexwithBS)                     ) ;
      store("muon_dxy"                             , muIt->globalTrack()->dxy(beamspot)                               ) ;
      store("muon_dsz"                             , muIt->globalTrack()->dsz(beamspot)                               ) ;
      store("muon_dz"                              , muIt->globalTrack()->dz(beamspot)                                ) ;
      store("muon_vx"                              , muIt->globalTrack()->vx()                                        ) ;
      store("muon_vy"                              , muIt->globalTrack()->vy()                                        ) ;
      store("muon_vz"                              , muIt->globalTrack()->vz()                                        ) ;
      store("muon_outerPt"                         , muIt->globalTrack()->outerPt()                                   ) ;
      store("muon_outerEta"                        , muIt->globalTrack()->outerEta()                                  ) ;
      store("muon_outerPhi"                        , muIt->globalTrack()->outerPhi()                                  ) ;
      store("muon_outerTheta"                      , muIt->globalTrack()->outerTheta()                                ) ;
      
      store("muon_qoverpError"                     , muIt->globalTrack()->qoverpError()                               ) ;
      store("muon_ptError"                         , muIt->globalTrack()->ptError()                                   ) ;
      store("muon_thetaError"                      , muIt->globalTrack()->thetaError()                                ) ;
      store("muon_lambdaError"                     , muIt->globalTrack()->lambdaError()                               ) ;
      store("muon_phiError"                        , muIt->globalTrack()->phiError()                                  ) ;
      store("muon_dxyError"                        , muIt->globalTrack()->dxyError()                                  ) ;
      store("muon_d0Error"                         , muIt->globalTrack()->d0Error()                                   ) ;
      store("muon_dszError"                        , muIt->globalTrack()->dszError()                                  ) ;
      store("muon_dzError"                         , muIt->globalTrack()->dzError()                                   ) ;
      
      store("muon_numberOfValidPixelHits"          , muIt->innerTrack()->hitPattern().numberOfValidPixelHits()        ) ;
      store("muon_numberOfValidTrackerHits"        , muIt->globalTrack()->hitPattern().numberOfValidTrackerHits()     ) ;                                  
      store("muon_numberOfValidMuonHits"           , muIt->globalTrack()->hitPattern().numberOfValidMuonHits()        ) ;                            
      store("muon_numberOfValidHits"               , muIt->globalTrack()->numberOfValidHits()                         ) ;        
      store("muon_trackerLayersWithMeasurement"    , muIt->track()->hitPattern().trackerLayersWithMeasurement()       ) ;        
      store("muon_numberOfLostHits"                , muIt->globalTrack()->numberOfLostHits()                          ) ;        
      store("muon_numberOfMatchedStations"         , muIt->numberOfMatchedStations()                                  ) ; 
      store("muon_validFraction"                   , muIt->globalTrack()->validFraction()                             ) ; 
             
      store("muon_isGlobalMuon"                    , muIt->isGlobalMuon()                                             ) ;        
      store("muon_isTrackerMuon"                   , muIt->isTrackerMuon()                                            ) ;        
      store("muon_isPFMuon"                        , muIt->isPFMuon()                                                 ) ;        
      store("muon_isPFIsolationValid"              , muIt->isPFIsolationValid()                                       ) ;        
      store("muon_chi2"                            , muIt->globalTrack()->chi2()                                      ) ;
      store("muon_ndof"                            , muIt->globalTrack()->ndof()                                      ) ;
      store("muon_normalizedChi2"                  , muIt->globalTrack()->normalizedChi2()                            ) ;
      store("muon_innerPosition_x"                 , muIt->globalTrack()->innerPosition().X()                         ) ;
      store("muon_innerPosition_y"                 , muIt->globalTrack()->innerPosition().Y()                         ) ;
      store("muon_innerPosition_z"                 , muIt->globalTrack()->innerPosition().Z()                         ) ;
      
      store("muon_isolationR03_sumPt"              , muIt->isolationR03().sumPt                                       ) ;
      store("muon_isolationR03_trackerVetoPt"      , muIt->isolationR03().trackerVetoPt                               ) ;
      store("muon_isolationR03_emEt"               , muIt->isolationR03().emEt                                        ) ;
      store("muon_isolationR03_emVetoEt"           , muIt->isolationR03().emVetoEt                                    ) ;
      store("muon_isolationR03_hadEt"              , muIt->isolationR03().hadEt                                       ) ;
      store("muon_isolationR03_hadVetoEt"          , muIt->isolationR03().hadVetoEt                                   ) ;
      
      store("muon_isolationR05_sumPt"              , muIt->isolationR05().sumPt                                       ) ;
      store("muon_isolationR05_trackerVetoPt"      , muIt->isolationR05().trackerVetoPt                               ) ;
      store("muon_isolationR05_emEt"               , muIt->isolationR05().emEt                                        ) ;
      store("muon_isolationR05_emVetoEt"           , muIt->isolationR05().emVetoEt                                    ) ;
      store("muon_isolationR05_hadEt"              , muIt->isolationR05().hadEt                                       ) ;
      store("muon_isolationR05_hadVetoEt"          , muIt->isolationR05().hadVetoEt                                   ) ;
      
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
      
      store("muon_tevOptimized_charge"             , tevOptimizedTrk.first->charge()                                  ) ;
      store("muon_tevOptimized_pt"                 , tevOptimizedTrk.first->pt()                                      ) ;
      store("muon_tevOptimized_eta"                , tevOptimizedTrk.first->eta()                                     ) ;
      store("muon_tevOptimized_phi"                , tevOptimizedTrk.first->phi()                                     ) ;
      store("muon_tevOptimized_theta"              , tevOptimizedTrk.first->theta()                                   ) ;
      store("muon_tevOptimized_px"                 , tevOptimizedTrk.first->px()                                      ) ;
      store("muon_tevOptimized_py"                 , tevOptimizedTrk.first->py()                                      ) ;
      store("muon_tevOptimized_pz"                 , tevOptimizedTrk.first->pz()                                      ) ;
      store("muon_tevOptimized_d0"                 , tevOptimizedTrk.first->d0()                                      ) ;
      store("muon_tevOptimized_dz"                 , tevOptimizedTrk.first->dz()                                      ) ;
      store("muon_tevOptimized_dz_beamSpot"        , tevOptimizedTrk.first->dz(beamspot)                              ) ;
      store("muon_tevOptimized_dz_firstPVtx"       , tevOptimizedTrk.first->dz(firstpvertex)                          ) ;
      store("muon_tevOptimized_dz_firstPVtxwithBS" , tevOptimizedTrk.first->dz(firstpvertexwithBS)                    ) ;
      store("muon_tevOptimized_dxy"                , tevOptimizedTrk.first->dxy()                                     ) ;
      store("muon_tevOptimized_dxy_beamSpot"       , tevOptimizedTrk.first->dxy(beamspot)                             ) ;
      store("muon_tevOptimized_dxy_firstPVtx"      , tevOptimizedTrk.first->dxy(firstpvertex)                         ) ;
      store("muon_tevOptimized_dxy_firstPVtxwithBS", tevOptimizedTrk.first->dxy(firstpvertexwithBS)                   ) ;
      
      store("muon_tevOptimized_ptError"            , tevOptimizedTrk.first->ptError()                                 ) ;
      store("muon_tevOptimized_etaError"           , tevOptimizedTrk.first->etaError()                                ) ;
      store("muon_tevOptimized_phiError"           , tevOptimizedTrk.first->phiError()                                ) ;
      store("muon_tevOptimized_thetaError"         , tevOptimizedTrk.first->thetaError()                              ) ;
      store("muon_tevOptimized_d0Error"            , tevOptimizedTrk.first->d0Error()                                 ) ;
      store("muon_tevOptimized_dzError"            , tevOptimizedTrk.first->dzError()                                 ) ;
      store("muon_tevOptimized_dxyError"           , tevOptimizedTrk.first->dxyError()                                ) ;

      std::string branchPrefix = "muMatch_" ;
      std::vector<std::string> filterNames ;
      filterNames.push_back("hltL1sMu16Eta2p1") ;
      filterNames.push_back("hltL1sL1Mu3p5EG12") ;
      filterNames.push_back("hltL1Mu3p5EG12L3Filtered22") ;
      filterNames.push_back("hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40Q") ;
      
      for(unsigned iFilter=0 ; iFilter<filterNames.size() ; iFilter++){
        bool muMatch = false ;
        // It is important to specify the right HLT process for the filter, not doing this is a common bug
        trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterNames.at(iFilter),"",trigEventTag.process())); 
        if(filterIndex<trigEvent->sizeFilters()){ 
          const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
          const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
          // Now loop over the trigger objects passing filter
          for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); ++keyIt) { 
            const trigger::TriggerObject& obj = trigObjColl[*keyIt];
            if(deltaR(tevOptimizedTrk.first->eta(), tevOptimizedTrk.first->phi(), obj.eta(), obj.phi())<1.){
              muMatch = true ;
            }
          }
        }//end filter size check
        std::string branchName = branchPrefix + filterNames.at(iFilter) ;
        store(branchName, muMatch) ;
      }
    }
  }
}

void IIHEModuleMuon::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMuon::beginEvent(){}
void IIHEModuleMuon::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMuon::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMuon);
