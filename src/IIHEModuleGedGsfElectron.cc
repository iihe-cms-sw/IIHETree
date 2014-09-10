#include "UserCode/IIHETree/interface/IIHEModuleGedGsfElectron.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

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

IIHEModuleGedGsfElectron::IIHEModuleGedGsfElectron(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleGedGsfElectron::~IIHEModuleGedGsfElectron(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleGedGsfElectron::beginJob(){
  addBranch("gsf_n", kUInt) ;
  addBranch("gsf_classification", kVectorInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("gsf_energy") ;
  addBranch("gsf_p") ;
  addBranch("gsf_pt") ;
  addBranch("gsf_scE1x5") ;
  addBranch("gsf_scE5x5") ;
  addBranch("gsf_scE2x5Max") ;
  addBranch("gsf_eta") ;
  addBranch("gsf_phi") ;
  addBranch("gsf_theta") ;
  addBranch("gsf_px") ;
  addBranch("gsf_py") ;
  addBranch("gsf_pz") ;
  addBranch("gsf_superClusterEta") ;
  addBranch("gsf_superClusterEnergy") ;
  addBranch("gsf_caloEnergy") ;
  addBranch("gsf_deltaEtaSuperClusterTrackAtVtx") ;
  addBranch("gsf_deltaPhiSuperClusterTrackAtVtx") ;
  addBranch("gsf_hadronicOverEm") ;
  addBranch("gsf_hcalDepth1OverEcal") ;
  addBranch("gsf_hcalDepth2OverEcal") ;
  addBranch("gsf_dr03TkSumPt") ;
  addBranch("gsf_dr03EcalRecHitSumEt") ;
  addBranch("gsf_dr03HcalDepth1TowerSumEt") ;
  addBranch("gsf_dr03HcalDepth2TowerSumEt") ;
  addBranch("gsf_charge", kVectorInt) ;
  addBranch("gsf_sigmaIetaIeta") ;
  addBranch("gsf_ecaldrivenSeed"   , kVectorBool) ;
  addBranch("gsf_trackerdrivenSeed", kVectorBool) ;
  addBranch("gsf_isEB", kVectorBool) ;
  addBranch("gsf_isEE", kVectorBool) ;
  addBranch("gsf_deltaEtaSeedClusterTrackAtCalo") ;
  addBranch("gsf_deltaPhiSeedClusterTrackAtCalo") ;
  addBranch("gsf_ecalEnergy") ;
  addBranch("gsf_eSuperClusterOverP") ;
  addBranch("gsf_dxy") ;
  addBranch("gsf_dxy_beamSpot") ;
  addBranch("gsf_dxy_firstPVtx") ;
  addBranch("gsf_dxy_firstPVtxwithBS") ;
  addBranch("gsf_dxyError") ;
  addBranch("gsf_dz") ;
  addBranch("gsf_dz_beamSpot") ;
  addBranch("gsf_dz_firstPVtx") ;
  addBranch("gsf_dz_firstPVtxwithBS") ;
  addBranch("gsf_dzError") ;
  addBranch("gsf_vz") ;
  addBranch("gsf_numberOfValidHits", kVectorInt) ;
  addBranch("gsf_nLostInnerHits"   , kVectorInt) ;
  addBranch("gsf_nLostOuterHits"   , kVectorInt) ;
  addBranch("gsf_convFlags"        , kVectorInt) ;
  addBranch("gsf_convDist") ;
  addBranch("gsf_convDcot") ;
  addBranch("gsf_convRadius") ;
  addBranch("gsf_fBrem") ;
  addBranch("gsf_e1x5") ;
  addBranch("gsf_e2x5Max") ;
  addBranch("gsf_e5x5") ;
  addBranch("gsf_hitsinfo", kVectorVectorInt) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleGedGsfElectron::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Delegate default electron collection name to IIHEAnalysis class
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  
  math::XYZPoint firstpvertexwithBS(0.0,0.0,0.0);
  Handle<reco::VertexCollection> primaryVertexCollwithBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVertexCollwithBS);
  const reco::VertexCollection* pvcollwithBS = primaryVertexCollwithBS.product();

  if(pvcollwithBS->size()>0){
    reco::VertexCollection::const_iterator firstpv = pvcollwithBS->begin();
    firstpvertexwithBS.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByLabel("offlineBeamSpot", theBeamSpot);
  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint firstpvertex(0.0,0.0,0.0);
  
  store("gsf_n", (unsigned int) electrons.size()) ;
  for(reco::GsfElectronCollection::const_iterator gsfiter=electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
    
    //Fill the gsf related variables
    int gsf_nLostInnerHits = gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
    int gsf_nLostOuterHits = gsfiter->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits() ;
    store("gsf_energy"                        , gsfiter->energy()) ;
    store("gsf_p"                             , gsfiter->p()) ;
    store("gsf_pt"                            , gsfiter->pt()) ;
    store("gsf_classification"                , gsfiter->classification()) ;
    store("gsf_scE1x5"                        , gsfiter->scE1x5()) ;
    store("gsf_scE5x5"                        , gsfiter->scE5x5()) ;
    store("gsf_scE2x5Max"                     , gsfiter->scE2x5Max()) ;
    store("gsf_eta"                           , gsfiter->eta()) ;
    store("gsf_phi"                           , gsfiter->phi()) ;
    store("gsf_theta"                         , gsfiter->theta());
    store("gsf_px"                            , gsfiter->px()) ;
    store("gsf_py"                            , gsfiter->py()) ;
    store("gsf_pz"                            , gsfiter->pz()) ;
    store("gsf_superClusterEta"               , gsfiter->superCluster()->eta()) ;
    store("gsf_superClusterEnergy"            , gsfiter->superCluster()->energy()) ;
    store("gsf_caloEnergy"                    , gsfiter->caloEnergy()) ;
    store("gsf_deltaEtaSuperClusterTrackAtVtx", gsfiter->deltaEtaSuperClusterTrackAtVtx()) ;
    store("gsf_deltaPhiSuperClusterTrackAtVtx", gsfiter->deltaPhiSuperClusterTrackAtVtx()) ;
    store("gsf_hadronicOverEm"                , gsfiter->hadronicOverEm()) ;
    store("gsf_hcalDepth1OverEcal"            , gsfiter->hcalDepth1OverEcal()) ;
    store("gsf_hcalDepth2OverEcal"            , gsfiter->hcalDepth2OverEcal()) ;
    store("gsf_dr03TkSumPt"                   , gsfiter->dr03TkSumPt()) ;
    store("gsf_dr03EcalRecHitSumEt"           , gsfiter->dr03EcalRecHitSumEt()) ;
    store("gsf_dr03HcalDepth1TowerSumEt"      , gsfiter->dr03HcalDepth1TowerSumEt()) ;
    store("gsf_dr03HcalDepth2TowerSumEt"      , gsfiter->dr03HcalDepth2TowerSumEt()) ;
    store("gsf_charge"                        , gsfiter->charge()) ;
    store("gsf_sigmaIetaIeta"                 , gsfiter->sigmaIetaIeta()) ;
    store("gsf_ecaldrivenSeed"                , gsfiter->ecalDrivenSeed()) ;
    store("gsf_trackerdrivenSeed"             , gsfiter->trackerDrivenSeed()) ;
    store("gsf_isEB"                          , gsfiter->isEB());
    store("gsf_isEE"                          , gsfiter->isEE());
    store("gsf_deltaEtaSeedClusterTrackAtCalo", gsfiter->deltaEtaSeedClusterTrackAtCalo());
    store("gsf_deltaPhiSeedClusterTrackAtCalo", gsfiter->deltaPhiSeedClusterTrackAtCalo());
    store("gsf_ecalEnergy"                    , gsfiter->ecalEnergy());
    store("gsf_eSuperClusterOverP"            , gsfiter->eSuperClusterOverP());
    store("gsf_dxy"                           , gsfiter->gsfTrack()->dxy());
    store("gsf_dxy_beamSpot"                  , gsfiter->gsfTrack()->dxy(beamspot));
    store("gsf_dxy_firstPVtx"                 , gsfiter->gsfTrack()->dxy(firstpvertex));
    store("gsf_dxy_firstPVtxwithBS"           , gsfiter->gsfTrack()->dxy(firstpvertexwithBS));
    store("gsf_dxyError"                      , gsfiter->gsfTrack()->dxyError());
    store("gsf_dz"                            , gsfiter->gsfTrack()->dz());
    store("gsf_dz_beamSpot"                   , gsfiter->gsfTrack()->dz(beamspot));
    store("gsf_dz_firstPVtx"                  , gsfiter->gsfTrack()->dz(firstpvertex));
    store("gsf_dz_firstPVtxwithBS"            , gsfiter->gsfTrack()->dz(firstpvertexwithBS));
    store("gsf_dzError"                       , gsfiter->gsfTrack()->dzError()); 
    store("gsf_vz"                            , gsfiter->gsfTrack()->vz());
    store("gsf_numberOfValidHits"             , gsfiter->gsfTrack()->numberOfValidHits());
    store("gsf_nLostInnerHits"                , gsf_nLostInnerHits);
    store("gsf_nLostOuterHits"                , gsf_nLostOuterHits);
    store("gsf_convFlags"                     , gsfiter->convFlags());
    store("gsf_convDist"                      , gsfiter->convDist());
    store("gsf_convDcot"                      , gsfiter->convDcot());
    store("gsf_convRadius"                    , gsfiter->convRadius());
    store("gsf_fBrem"                         , gsfiter->fbrem());
    store("gsf_e1x5"                          , gsfiter->e1x5()) ;
    store("gsf_e2x5Max"                       , gsfiter->e2x5Max()) ;
    store("gsf_e5x5"                          , gsfiter->e5x5()) ;
    
    //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/TrackReco/interface/HitPattern.h?revision=1.32&view=markup
    reco::HitPattern kfHitPattern = gsfiter->gsfTrack()->hitPattern();
    int nbtrackhits = kfHitPattern.numberOfHits();
    std::vector<int> gsf_hitsinfo ;
    for(int hititer=0 ; hititer<25 ; hititer++){
      int myhitbin = (hititer<nbtrackhits) ? kfHitPattern.getHitPattern(hititer) : 0 ;
      gsf_hitsinfo.push_back(myhitbin) ;
    }
    store("gsf_hitsinfo", gsf_hitsinfo) ;
  }
}

void IIHEModuleGedGsfElectron::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleGedGsfElectron::beginEvent(){}
void IIHEModuleGedGsfElectron::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleGedGsfElectron::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleGedGsfElectron);
