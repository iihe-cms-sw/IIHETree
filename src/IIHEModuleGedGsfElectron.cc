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

IIHEModuleGedGsfElectron::IIHEModuleGedGsfElectron(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
  ETThreshold_ = iConfig.getUntrackedParameter<double>("electrons_ETThreshold", 0.0 ) ;
}
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
  setBranchType(kVectorBool) ;
  addBranch("gsf_ecaldrivenSeed"   ) ;
  addBranch("gsf_trackerdrivenSeed") ;
  addBranch("gsf_isEB") ;
  addBranch("gsf_isEE") ;
  setBranchType(kVectorFloat) ;
  addBranch("gsf_deltaEtaSeedClusterTrackAtCalo") ;
  addBranch("gsf_deltaPhiSeedClusterTrackAtCalo") ;
  addBranch("gsf_ecalEnergy") ;
  addBranch("gsf_eSuperClusterOverP") ;
  addBranch("gsf_dxy") ;
  addBranch("gsf_dxy_beamSpot") ;
  addBranch("gsf_dxy_firstPVtx") ;
  addBranch("gsf_dxyError") ;
  addBranch("gsf_dz") ;
  addBranch("gsf_dz_beamSpot") ;
  addBranch("gsf_dz_firstPVtx") ;
  addBranch("gsf_dzError") ;
  addBranch("gsf_vz") ;
  setBranchType(kVectorInt) ;
  addBranch("gsf_numberOfValidHits") ;
  addBranch("gsf_nLostInnerHits"   ) ;
  addBranch("gsf_nLostOuterHits"   ) ;
  addBranch("gsf_convFlags"        ) ;
  setBranchType(kVectorFloat) ;
  addBranch("gsf_convDist") ;
  addBranch("gsf_convDcot") ;
  addBranch("gsf_convRadius") ;
  addBranch("gsf_fBrem") ;
  addBranch("gsf_e1x5") ;
  addBranch("gsf_e2x5Max") ;
  addBranch("gsf_e5x5") ;
  addBranch("gsf_r9") ;
  addBranch("gsf_hitsinfo", kVectorVectorInt) ;
  
//CHOOSE_RELEASE_START CMSSW_7_4_4
  setBranchType(kVectorFloat) ;
  addBranch("gsf_pixelMatch_dPhi1") ;
  addBranch("gsf_pixelMatch_dPhi2") ;
  addBranch("gsf_pixelMatch_dRz1" ) ;
  addBranch("gsf_pixelMatch_dRz2" ) ;
  setBranchType(kVectorInt) ;
  addBranch("gsf_pixelMatch_subDetector1") ;
  addBranch("gsf_pixelMatch_subDetector2") ;
//CHOOSE_RELEASE_END CMSSW_7_4_4
/*CHOOSE_RELEASE_START CMSSW_7_2_0 CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11
CHOOSE_RELEASE_END CMSSW_7_2_0 CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11*/
  
  addBranch("gsf_mc_bestDR", kVectorFloat) ;
  addBranch("gsf_mc_index" , kVectorInt  ) ;
  addBranch("gsf_mc_ERatio", kVectorFloat) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleGedGsfElectron::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Delegate default electron collection name to IIHEAnalysis class
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  
  math::XYZPoint* beamspot     = parent_->getBeamspot() ;
  math::XYZPoint* firstpvertex = parent_->getFirstPrimaryVertex() ;
  
  unsigned int gsf_n = 0 ;
  for(reco::GsfElectronCollection::const_iterator gsfiter=electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
    
    float ET = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
    if(ET<ETThreshold_ && gsfiter->pt()<ETThreshold_) continue ;
    gsf_n++ ;
    
    //Fill the gsf related variables

//CHOOSE_RELEASE_START DEFAULT CMSSW_7_4_4 CMSSW_7_3_0 CMSSW_7_2_0
    int gsf_nLostInnerHits = gsfiter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) ;
    int gsf_nLostOuterHits = gsfiter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS) ;
//CHOOSE_RELEASE_END DEFAULT CMSSW_7_4_4 CMSSW_7_3_0 CMSSW_7_2_0
/*CHOOSE_RELEASE_START CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11
    int gsf_nLostInnerHits = gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
    int gsf_nLostOuterHits = gsfiter->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits() ;
CHOOSE_RELEASE_END CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11*/
    
    store("gsf_energy"                        , gsfiter->energy()                        ) ;
    store("gsf_p"                             , gsfiter->p()                             ) ;
    store("gsf_pt"                            , gsfiter->pt()                            ) ;
    store("gsf_classification"                , gsfiter->classification()                ) ;
    store("gsf_scE1x5"                        , gsfiter->scE1x5()                        ) ;
    store("gsf_scE5x5"                        , gsfiter->scE5x5()                        ) ;
    store("gsf_scE2x5Max"                     , gsfiter->scE2x5Max()                     ) ;
    store("gsf_eta"                           , gsfiter->eta()                           ) ;
    store("gsf_phi"                           , gsfiter->phi()                           ) ;
    store("gsf_theta"                         , gsfiter->theta()                         ) ;
    store("gsf_px"                            , gsfiter->px()                            ) ;
    store("gsf_py"                            , gsfiter->py()                            ) ;
    store("gsf_pz"                            , gsfiter->pz()                            ) ;
    store("gsf_superClusterEta"               , gsfiter->superCluster()->eta()           ) ;
    store("gsf_superClusterEnergy"            , gsfiter->superCluster()->energy()        ) ;
    store("gsf_caloEnergy"                    , gsfiter->caloEnergy()                    ) ;
    store("gsf_deltaEtaSuperClusterTrackAtVtx", gsfiter->deltaEtaSuperClusterTrackAtVtx()) ;
    store("gsf_deltaPhiSuperClusterTrackAtVtx", gsfiter->deltaPhiSuperClusterTrackAtVtx()) ;
    store("gsf_hadronicOverEm"                , gsfiter->hadronicOverEm()                ) ;
    store("gsf_hcalDepth1OverEcal"            , gsfiter->hcalDepth1OverEcal()            ) ;
    store("gsf_hcalDepth2OverEcal"            , gsfiter->hcalDepth2OverEcal()            ) ;
    store("gsf_dr03TkSumPt"                   , gsfiter->dr03TkSumPt()                   ) ;
    store("gsf_dr03EcalRecHitSumEt"           , gsfiter->dr03EcalRecHitSumEt()           ) ;
    store("gsf_dr03HcalDepth1TowerSumEt"      , gsfiter->dr03HcalDepth1TowerSumEt()      ) ;
    store("gsf_dr03HcalDepth2TowerSumEt"      , gsfiter->dr03HcalDepth2TowerSumEt()      ) ;
    store("gsf_charge"                        , gsfiter->charge()                        ) ;
    store("gsf_sigmaIetaIeta"                 , gsfiter->sigmaIetaIeta()                 ) ;
    store("gsf_ecaldrivenSeed"                , gsfiter->ecalDrivenSeed()                ) ;
    store("gsf_trackerdrivenSeed"             , gsfiter->trackerDrivenSeed()             ) ;
    store("gsf_isEB"                          , gsfiter->isEB()                          ) ;
    store("gsf_isEE"                          , gsfiter->isEE()                          ) ;
    store("gsf_deltaEtaSeedClusterTrackAtCalo", gsfiter->deltaEtaSeedClusterTrackAtCalo()) ;
    store("gsf_deltaPhiSeedClusterTrackAtCalo", gsfiter->deltaPhiSeedClusterTrackAtCalo()) ;
    store("gsf_ecalEnergy"                    , gsfiter->ecalEnergy()                    ) ;
    store("gsf_eSuperClusterOverP"            , gsfiter->eSuperClusterOverP()            ) ;
    store("gsf_dxy"                           , gsfiter->gsfTrack()->dxy()               ) ;
    store("gsf_dxy_beamSpot"                  , gsfiter->gsfTrack()->dxy(*beamspot)      ) ;
    store("gsf_dxy_firstPVtx"                 , gsfiter->gsfTrack()->dxy(*firstpvertex)  ) ;
    store("gsf_dxyError"                      , gsfiter->gsfTrack()->dxyError()          ) ;
    store("gsf_dz"                            , gsfiter->gsfTrack()->dz()                ) ;
    store("gsf_dz_beamSpot"                   , gsfiter->gsfTrack()->dz(*beamspot)       ) ;
    store("gsf_dz_firstPVtx"                  , gsfiter->gsfTrack()->dz(*firstpvertex)   ) ;
    store("gsf_dzError"                       , gsfiter->gsfTrack()->dzError()           ) ; 
    store("gsf_vz"                            , gsfiter->gsfTrack()->vz()                ) ;
    store("gsf_numberOfValidHits"             , gsfiter->gsfTrack()->numberOfValidHits() ) ;
    store("gsf_nLostInnerHits"                , gsf_nLostInnerHits                       ) ;
    store("gsf_nLostOuterHits"                , gsf_nLostOuterHits                       ) ;
    store("gsf_convFlags"                     , gsfiter->convFlags()                     ) ;
    store("gsf_convDist"                      , gsfiter->convDist()                      ) ;
    store("gsf_convDcot"                      , gsfiter->convDcot()                      ) ;
    store("gsf_convRadius"                    , gsfiter->convRadius()                    ) ;
    store("gsf_fBrem"                         , gsfiter->fbrem()                         ) ;
    store("gsf_e1x5"                          , gsfiter->e1x5()                          ) ;
    store("gsf_r9"                            , gsfiter->r9()                            ) ;
    store("gsf_e2x5Max"                       , gsfiter->e2x5Max()                       ) ;
    store("gsf_e5x5"                          , gsfiter->e5x5()                          ) ;
    
    //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/TrackReco/interface/HitPattern.h?revision=1.32&view=markup
    reco::HitPattern kfHitPattern = gsfiter->gsfTrack()->hitPattern();

//CHOOSE_RELEASE_START DEFAULT CMSSW_7_4_4 CMSSW_7_3_0 CMSSW_7_2_0
    int nbtrackhits = kfHitPattern.numberOfHits(reco::HitPattern::TRACK_HITS) ;
//CHOOSE_RELEASE_END DEFAULT CMSSW_7_4_4 CMSSW_7_3_0 CMSSW_7_2_0
/*CHOOSE_RELEASE_START CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11
    int nbtrackhits = kfHitPattern.numberOfHits() ;
CHOOSE_RELEASE_END CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11*/

    std::vector<int> gsf_hitsinfo ;
    for(int hititer=0 ; hititer<25 ; hititer++){
      
//CHOOSE_RELEASE_START CMSSW_7_4_4 CMSSW_7_3_0 CMSSW_7_2_0
      int myhitbin = (hititer<nbtrackhits) ? kfHitPattern.getHitPattern(reco::HitPattern::TRACK_HITS, hititer) : 0 ;
//CHOOSE_RELEASE_END CMSSW_7_4_4 CMSSW_7_3_0 CMSSW_7_2_0
/*CHOOSE_RELEASE_START CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11
      int myhitbin = (hititer<nbtrackhits) ? kfHitPattern.getHitPattern(hititer) : 0 ;
CHOOSE_RELEASE_END CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11*/
      
      gsf_hitsinfo.push_back(myhitbin) ;
    }
    store("gsf_hitsinfo", gsf_hitsinfo) ;
    
//CHOOSE_RELEASE_START CMSSW_7_4_4
    store("gsf_pixelMatch_dPhi1"       , gsfiter->pixelMatchDPhi1()       ) ;
    store("gsf_pixelMatch_dPhi2"       , gsfiter->pixelMatchDPhi2()       ) ;
    store("gsf_pixelMatch_dRz1"        , gsfiter->pixelMatchDRz1()        ) ;
    store("gsf_pixelMatch_dRz2"        , gsfiter->pixelMatchDRz2()        ) ;
    store("gsf_pixelMatch_subDetector1", gsfiter->pixelMatchSubdetector1()) ;
    store("gsf_pixelMatch_subDetector2", gsfiter->pixelMatchSubdetector2()) ;
//CHOOSE_RELEASE_END CMSSW_7_4_4
/*CHOOSE_RELEASE_START CMSSW_7_2_0 CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11
CHOOSE_RELEASE_END CMSSW_7_2_0 CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11*/
    
    // Now apply truth matching.
    int index = MCTruth_matchEtaPhi_getIndex(gsfiter->eta(), gsfiter->phi()) ;
    if(index>=0){
      const MCTruthObject* MCTruth = MCTruth_getRecordByIndex(index) ;
      store("gsf_mc_bestDR", deltaR(gsfiter->eta(), gsfiter->phi(), MCTruth->eta(), MCTruth->phi())) ;
      store("gsf_mc_index" , index) ;
      store("gsf_mc_ERatio", gsfiter->energy()/MCTruth->energy()) ;
    }
    else{
      store("gsf_mc_bestDR", 999.0) ;
      store("gsf_mc_index" ,    -1) ;
      store("gsf_mc_ERatio", 999.0) ;
    }
  }
  store("gsf_n", gsf_n) ;
}

void IIHEModuleGedGsfElectron::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleGedGsfElectron::beginEvent(){}
void IIHEModuleGedGsfElectron::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleGedGsfElectron::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleGedGsfElectron);
