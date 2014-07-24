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
  addBranch("muon_n", kInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("muon_pt") ;
  addBranch("muon_ptError") ;
  addBranch("muon_gTrk_pt") ;
  addBranch("muon_gTrk_ptError") ;
  addBranch("muon_eta") ;
  addBranch("muon_etaError") ;
  addBranch("muon_phi") ;
  addBranch("muon_phiError") ;
  addBranch("muon_theta") ;
  addBranch("muon_thetaError") ;
  addBranch("muon_outerPt") ;
  addBranch("muon_outerEta") ;
  addBranch("muon_outerPhi") ;
  addBranch("muon_outerTheta") ;
  addBranch("muon_px") ;
  addBranch("muon_py") ;
  addBranch("muon_pz") ;
  setBranchType(kVectorInt) ;
  addBranch("muon_charge") ;
  addBranch("muon_nhitspixel") ;
  addBranch("muon_nhitstrack") ;
  addBranch("muon_nhitsmuons") ;
  addBranch("muon_nhitstotal") ;
  addBranch("muon_nlayerswithhits") ;
  addBranch("muon_nlosthits") ;
  addBranch("muon_nSegmentMatch") ;
  setBranchType(kVectorBool) ;
  addBranch("muon_isTrackerMuon") ;
  addBranch("muon_isPFMuon") ;
  addBranch("muon_isPFIsolationValid") ;
  addBranch("muon_chi2"    , kVectorFloat) ;
  addBranch("muon_ndof"    , kVectorInt  ) ;
  addBranch("muon_normChi2", kVectorFloat) ;
  setBranchType(kVectorFloat) ;
  addBranch("muon_d0") ;
  addBranch("muon_d0Error") ;
  addBranch("muon_dz_cmsCenter") ;
  addBranch("muon_dz_beamSpot") ;
  addBranch("muon_dz_firstPVtx") ;
  addBranch("muon_dz_firstPVtxwithBS") ;
  addBranch("muon_dzError") ;
  addBranch("muon_dxy_cmsCenter") ;
  addBranch("muon_dxy_beamSpot") ;
  addBranch("muon_dxy_firstPVtx") ;
  addBranch("muon_dxy_firstPVtxwithBS") ;
  addBranch("muon_dxyError") ;
  addBranch("muon_trackIso03") ;
  addBranch("muon_trackIso05") ;
  addBranch("muon_trackIso03_ptInVeto") ;
  addBranch("muon_trackIso05_ptInVeto") ;
  addBranch("muon_emIso03") ;
  addBranch("muon_emIso05") ;
  addBranch("muon_emIso03_ptInVeto") ;
  addBranch("muon_emIso05_ptInVeto") ;
  addBranch("muon_hadIso03") ;
  addBranch("muon_hadIso05") ;
  addBranch("muon_hadIso03_ptInVeto") ;
  addBranch("muon_hadIso05_ptInVeto") ;
  addBranch("muon_innerPosx") ;
  addBranch("muon_innerPosy") ;
  addBranch("muon_innerPosz") ;
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
  edm::Handle<reco::MuonCollection> muonCollection;
  iEvent.getByLabel("muons",muonCollection);
  const reco::MuonCollection* muons = muonCollection.product();
  store("muon_n", (int)muons->size()) ;
  for(reco::MuonCollection::const_iterator muIt = muons->begin(); muIt != muons->end(); ++muIt) {
    if (muIt->isGlobalMuon()) {
      // get TeV optimized track
      reco::Muon::MuonTrackTypePair tevOptMuTrk = muon::tevOptimized(*muIt, 200, 17., 40., 0.25) ;
      store("muon_pt"                 , tevOptMuTrk.first->pt()                                      ) ;
      store("muon_ptError"            , tevOptMuTrk.first->ptError()                                 ) ;
      store("muon_gTrk_pt"            , muIt->globalTrack()->pt()                                    ) ;
      store("muon_gTrk_ptError"       , muIt->globalTrack()->ptError()                               ) ;
      store("muon_eta"                , tevOptMuTrk.first->eta()                                     ) ;
      store("muon_etaError"           , tevOptMuTrk.first->etaError()                                ) ;
      store("muon_phi"                , tevOptMuTrk.first->phi()                                     ) ;
      store("muon_phiError"           , tevOptMuTrk.first->phiError()                                ) ;
      store("muon_theta"              , tevOptMuTrk.first->theta()                                   ) ;
      store("muon_thetaError"         , tevOptMuTrk.first->thetaError()                              ) ;
      store("muon_outerPt"            , muIt->globalTrack()->outerPt()                               ) ;
      store("muon_outerEta"           , muIt->globalTrack()->outerEta()                              ) ;
      store("muon_outerPhi"           , muIt->globalTrack()->outerPhi()                              ) ;
      store("muon_outerTheta"         , muIt->globalTrack()->outerTheta()                            ) ;
      store("muon_px"                 , tevOptMuTrk.first->px()                                      ) ;
      store("muon_py"                 , tevOptMuTrk.first->py()                                      ) ;
      store("muon_pz"                 , tevOptMuTrk.first->pz()                                      ) ;
      store("muon_charge"             , tevOptMuTrk.first->charge()                                  ) ;
      store("muon_nhitspixel"         , muIt->innerTrack()->hitPattern().numberOfValidPixelHits()    ) ;
      store("muon_nhitstrack"         , muIt->globalTrack()->hitPattern().numberOfValidTrackerHits() ) ;                                  
      store("muon_nhitsmuons"         , muIt->globalTrack()->hitPattern().numberOfValidMuonHits()    ) ;                            
      store("muon_nhitstotal"         , muIt->globalTrack()->numberOfValidHits()                     ) ;        
      store("muon_nlayerswithhits"    , muIt->track()->hitPattern().trackerLayersWithMeasurement()   ) ;        
      store("muon_nlosthits"          , muIt->globalTrack()->numberOfLostHits()                      ) ;        
      store("muon_nSegmentMatch"      , muIt->numberOfMatchedStations()                              ) ;        
      store("muon_isTrackerMuon"      , muIt->isTrackerMuon()                                        ) ;        
      store("muon_isPFMuon"           , muIt->isPFMuon()                                             ) ;        
      store("muon_isPFIsolationValid" , muIt->isPFIsolationValid()                                   ) ;        
      store("muon_chi2"               , muIt->globalTrack()->chi2()                                  ) ;
      store("muon_ndof"               , (int)(muIt->globalTrack()->ndof())                           ) ;
      store("muon_normChi2"           , muIt->globalTrack()->normalizedChi2()                        ) ;        
      store("muon_d0"                 , tevOptMuTrk.first->d0()                                      ) ;
      store("muon_d0Error"            , tevOptMuTrk.first->d0Error()                                 ) ;
      store("muon_dz_cmsCenter"       , tevOptMuTrk.first->dz()                                      ) ;
      store("muon_dz_beamSpot"        , tevOptMuTrk.first->dz(beamspot)                              ) ;
      store("muon_dz_firstPVtx"       , tevOptMuTrk.first->dz(firstpvertex)                          ) ;
      store("muon_dz_firstPVtxwithBS" , tevOptMuTrk.first->dz(firstpvertexwithBS)                    ) ;
      store("muon_dzError"            , tevOptMuTrk.first->dzError()                                 ) ;
      store("muon_dxy_cmsCenter"      , tevOptMuTrk.first->dxy()                                     ) ;
      store("muon_dxy_beamSpot"       , tevOptMuTrk.first->dxy(beamspot)                             ) ;
      store("muon_dxy_firstPVtx"      , tevOptMuTrk.first->dxy(firstpvertex)                         ) ;
      store("muon_dxy_firstPVtxwithBS", tevOptMuTrk.first->dxy(firstpvertexwithBS)                   ) ;
      store("muon_dxyError"           , tevOptMuTrk.first->dxyError()                                ) ;
      store("muon_innerPosx"          , muIt->globalTrack()->innerPosition().X()                     ) ;
      store("muon_innerPosy"          , muIt->globalTrack()->innerPosition().Y()                     ) ;
      store("muon_innerPosz"          , muIt->globalTrack()->innerPosition().Z()                     ) ;
      store("muon_trackIso03"         , muIt->isolationR03().sumPt                                   ) ;
      store("muon_trackIso05"         , muIt->isolationR05().sumPt                                   ) ;
      store("muon_trackIso03_ptInVeto", muIt->isolationR03().trackerVetoPt                           ) ;
      store("muon_trackIso05_ptInVeto", muIt->isolationR05().trackerVetoPt                           ) ;
      store("muon_emIso03"            , muIt->isolationR03().emEt                                    ) ;
      store("muon_emIso05"            , muIt->isolationR05().emEt                                    ) ;
      store("muon_emIso03_ptInVeto"   , muIt->isolationR03().emVetoEt                                ) ;
      store("muon_emIso05_ptInVeto"   , muIt->isolationR05().emVetoEt                                ) ;
      store("muon_hadIso03"           , muIt->isolationR03().hadEt                                   ) ;
      store("muon_hadIso05"           , muIt->isolationR05().hadEt                                   ) ;
      store("muon_hadIso03_ptInVeto"  , muIt->isolationR03().hadVetoEt                               ) ;
      store("muon_hadIso05_ptInVeto"  , muIt->isolationR05().hadVetoEt                               ) ;

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
            if(deltaR(tevOptMuTrk.first->eta(), tevOptMuTrk.first->phi(), obj.eta(), obj.phi())<1.){
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
