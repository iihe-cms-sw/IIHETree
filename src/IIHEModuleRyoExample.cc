#include "UserCode/IIHETree/interface/IIHEModuleRyoExample.h"

#include <TLorentzVector.h>

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleRyoExample::IIHEModuleRyoExample(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleRyoExample::~IIHEModuleRyoExample(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleRyoExample::beginJob(){
  addBranch("ryo_muon_n", kUInt) ;
  setBranchType(kVectorBool) ;
  addBranch("ryo_muonpass_Pt") ;
  addBranch("ryo_muonpass_Eta") ;
  addBranch("ryo_muonpass_IsGlobal") ;
  addBranch("ryo_muonpass_isTracker") ;
  addBranch("ryo_muonpass_normalizedChi2") ;
  addBranch("ryo_muonpass_numberOfValidMuonHits") ;
  addBranch("ryo_muonpass_numberOfMatchedStations") ;
  addBranch("ryo_muonpass_dxy") ;
  addBranch("ryo_muonpass_dz") ;
  addBranch("ryo_muonpass_numberOfValidPixelHits") ;
  addBranch("ryo_muonpass_trackerLayerWithMeasurement") ;
  addBranch("ryo_muonpass_isolation") ;
  addBranch("ryo_muonpass_all") ;
  
  addBranch("ryo_Zprime_n", kUInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("ryo_Zprime_pt") ;
  addBranch("ryo_Zprime_eta") ;
  addBranch("ryo_Zprime_phi") ;
  addBranch("ryo_Zprime_E") ;
  addBranch("ryo_Zprime_m") ;
  setBranchType(kVectorInt) ;
  addBranch("ryo_Zprime_mu1index") ;
  addBranch("ryo_Zprime_mu2index") ;
}

// ------------ method called to for each event  ------------
void IIHEModuleRyoExample::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Muon collections
  edm::Handle<reco::MuonCollection> muonCollection;
  iEvent.getByLabel("muons",muonCollection);
  const reco::MuonCollection* muons = muonCollection.product();
  unsigned int nMuon = 0 ;
  vector<pair<reco::Muon,int> > good_muons ;
  int muCounter = 0 ;
  for(reco::MuonCollection::const_iterator muIt = muons->begin(); muIt != muons->end(); ++muIt){
    bool muonpass_IsGlobal                    = muIt->isGlobalMuon() ;
    bool muonpass_isTracker                   = muIt->isTrackerMuon() ;
    bool muonpass_numberOfMatchedStations     = (muIt->numberOfMatchedStations() > 1) ;
    bool muonpass_pt                          = false ;
    bool muonpass_eta                         = false ;
    bool muonpass_normalizedChi2              = false ;
    bool muonpass_numberOfValidMuonHits       = false ;
    bool muonpass_dxy                         = false ;
    bool muonpass_dz                          = false ;
    bool muonpass_numberOfValidPixelHits      = false ;
    bool muonpass_trackerLayerWithMeasurement = false ;
    bool muonpass_isolation                   = false ;
    if(muonpass_IsGlobal){ // We need to protect against the globalTrack not existing for some cuts
      muonpass_pt                             = (muIt->globalTrack()->pt()>7) ;
      muonpass_eta                            = (fabs(muIt->globalTrack()->eta())<2.4) ;
      muonpass_normalizedChi2                 = (muIt->globalTrack()->normalizedChi2()<10) ;
      muonpass_numberOfValidMuonHits          = (muIt->globalTrack()->hitPattern().numberOfValidMuonHits() > 0) ;
      muonpass_isolation                      = ((muIt->isolationR03().sumPt/muIt->globalTrack()->pt()) < 0.10) ;
    }
    if(muonpass_isTracker){ // We need to protect against the globalTrack not existing for some cuts
      muonpass_dxy                            = (fabs(muIt->muonBestTrack()->dxy()) < 0.2) ;
      muonpass_dz                             = (fabs(muIt->muonBestTrack()->dz())  < 0.5) ;
      muonpass_numberOfValidPixelHits         = (muIt->innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ;
      muonpass_trackerLayerWithMeasurement    = (muIt->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 ) ;
    }
    bool muonpass_all                         = (muonpass_pt && muonpass_eta && muonpass_IsGlobal && muonpass_isTracker && muonpass_normalizedChi2 && muonpass_numberOfValidMuonHits && muonpass_numberOfMatchedStations && muonpass_dxy && muonpass_dz && muonpass_numberOfValidPixelHits && muonpass_trackerLayerWithMeasurement && muonpass_isolation) ;
    
    store("ryo_muonpass_Pt"                         , muonpass_pt                         ) ;
    store("ryo_muonpass_Eta"                        , muonpass_eta                        ) ;
    store("ryo_muonpass_IsGlobal"                   , muonpass_IsGlobal                   ) ;
    store("ryo_muonpass_isTracker"                  , muonpass_isTracker                  ) ;
    store("ryo_muonpass_normalizedChi2"             , muonpass_normalizedChi2             ) ;
    store("ryo_muonpass_numberOfValidMuonHits"      , muonpass_numberOfValidMuonHits      ) ;
    store("ryo_muonpass_numberOfMatchedStations"    , muonpass_numberOfMatchedStations    ) ;
    store("ryo_muonpass_dxy"                        , muonpass_dxy                        ) ;
    store("ryo_muonpass_dz"                         , muonpass_dz                         ) ;
    store("ryo_muonpass_numberOfValidPixelHits"     , muonpass_numberOfValidPixelHits     ) ;
    store("ryo_muonpass_trackerLayerWithMeasurement", muonpass_trackerLayerWithMeasurement) ;
    store("ryo_muonpass_isolation"                  , muonpass_isolation                  ) ;
    store("ryo_muonpass_all"                        , muonpass_all                        ) ;
    if(muonpass_all){
      nMuon++ ;
      good_muons.push_back(pair<reco::Muon,int>(*muIt,muCounter)) ;
    }
    muCounter++ ;
  }
  
  // Now take the good muons and make Zprime bosons from them
  // There are probably more dignified ways of doing it than this
  unsigned int nZprime = 0 ;
  for(unsigned int i=0 ; i<good_muons.size() ; i++){
    TLorentzVector m1p4 ;
    Muon mu1 = good_muons.at(i).first ;
    float pt1  = mu1.globalTrack()->pt() ;
    float eta1 = mu1.globalTrack()->eta() ;
    float phi1 = mu1.globalTrack()->phi() ;
    float m1   = 0.105 ;
    m1p4.SetPtEtaPhiM(pt1,eta1,phi1,m1) ;
    for(unsigned int j=i+1 ; j<good_muons.size() ; j++){
      TLorentzVector m2p4 ;
      Muon mu2 = good_muons.at(j).first ;
      float pt2  = mu2.globalTrack()->pt() ;
      float eta2 = mu2.globalTrack()->eta() ;
      float phi2 = mu2.globalTrack()->phi() ;
      float m2   = 0.105 ;
      m2p4.SetPtEtaPhiM(pt2,eta2,phi2,m2) ;
      TLorentzVector Zp4 = m1p4+m2p4 ;
      store("ryo_Zprime_pt" , Zp4.Pt() ) ;
      store("ryo_Zprime_eta", Zp4.Eta()) ;
      store("ryo_Zprime_phi", Zp4.Phi()) ;
      store("ryo_Zprime_E"  , Zp4.E()  ) ;
      store("ryo_Zprime_m"  , Zp4.M()  ) ;
      store("ryo_Zprime_mu1index", good_muons.at(i).second) ;
      store("ryo_Zprime_mu2index", good_muons.at(j).second) ;
      nZprime++ ;
    }
  }
  store("ryo_muon_n"  , nMuon  ) ;
  store("ryo_Zprime_n", nZprime) ;
}

void IIHEModuleRyoExample::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleRyoExample::beginEvent(){}
void IIHEModuleRyoExample::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleRyoExample::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleRyoExample);
