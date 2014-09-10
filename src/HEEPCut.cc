#include "UserCode/IIHETree/interface/HEEPCut.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                      Base class                                      //
//////////////////////////////////////////////////////////////////////////////////////////
HEEPCutBase::HEEPCutBase(std::string name){
  name_ = name ;
  index_ = -1 ;
  reset() ;
}
HEEPCutBase::~HEEPCutBase(){}
void HEEPCutBase::setStatus(bool value){ status_ = value ; }
bool HEEPCutBase::getStatus(){ return status_ ; }
void HEEPCutBase::reset(){ status_ = true ; }
bool HEEPCutBase::addBranch(IIHEModuleHEEP* mod){ return mod->addBranch(name_, kVectorBool) ; }
void HEEPCutBase::store(IIHEModuleHEEP* mod){ mod->store(name_, status_) ; }
bool HEEPCutBase::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){ return status_ ; }
bool HEEPCutBase::isBarrel(reco::GsfElectronCollection::const_iterator gsfiter){
  return (fabs(gsfiter->superCluster()->eta()) < 1.442) ;
}
bool HEEPCutBase::isEndcap(reco::GsfElectronCollection::const_iterator gsfiter){
  float eta = gsfiter->superCluster()->eta() ;
  return (fabs(eta) > 1.56 && fabs(eta) < 2.5) ;
}
int  HEEPCutBase::detectorRegion(reco::GsfElectronCollection::const_iterator gsfiter){
  if(isBarrel(gsfiter)) return kBarrel ;
  if(isEndcap(gsfiter)) return kEndcap ;
  return kNone ;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                     HEEP 4.1 cuts                                    //
//////////////////////////////////////////////////////////////////////////////////////////
HEEPCut_41_Et              ::HEEPCut_41_Et             (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_pt              ::HEEPCut_41_pt             (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_dEta            ::HEEPCut_41_dEta           (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_crack           ::HEEPCut_41_crack          (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_dEtaIn          ::HEEPCut_41_dEtaIn         (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_dPhiIn          ::HEEPCut_41_dPhiIn         (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_HOverE          ::HEEPCut_41_HOverE         (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_SigmaIetaIeta   ::HEEPCut_41_SigmaIetaIeta  (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_E2x5OverE5x5    ::HEEPCut_41_E2x5OverE5x5   (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_isolEMHadDepth1 ::HEEPCut_41_isolEMHadDepth1(std::string name): HEEPCutBase(name){} ;
HEEPCut_41_HcalIso2        ::HEEPCut_41_HcalIso2       (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_IsolPtTrks      ::HEEPCut_41_IsolPtTrks     (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_EcalDriven      ::HEEPCut_41_EcalDriven     (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_invalid         ::HEEPCut_41_invalid        (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_missingHits     ::HEEPCut_41_missingHits    (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_conversion      ::HEEPCut_41_conversion     (std::string name): HEEPCutBase(name){} ;
HEEPCut_41_dxyFirstPV      ::HEEPCut_41_dxyFirstPV     (std::string name): HEEPCutBase(name){} ;

bool HEEPCut_41_Et             ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  float Et  = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (Et > 35.0) ; break ;
    case kEndcap: result = (Et > 35.0) ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_pt             ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  bool result = true ;
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_dEta           ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  bool result = true ;
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_crack          ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  bool result = true ;
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_dEtaIn         ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  float deltaEta = gsfiter->deltaEtaSuperClusterTrackAtVtx() ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (fabs(deltaEta) < 0.005) ; break ;
    case kEndcap: result = (fabs(deltaEta) < 0.007) ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_dPhiIn         ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  float deltaPhi = gsfiter->deltaPhiSuperClusterTrackAtVtx() ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (fabs(deltaPhi) < 0.06) ; break ;
    case kEndcap: result = (fabs(deltaPhi) < 0.06) ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_HOverE         ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  float HOverE = gsfiter->hadronicOverEm() ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (HOverE < 0.05) ; break ;
    case kEndcap: result = (HOverE < 0.05) ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_SigmaIetaIeta  ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = true ; break ;
    case kEndcap: result = (gsfiter->sigmaIetaIeta() < 0.03) ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_E2x5OverE5x5   ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (gsfiter->scE2x5Max()/gsfiter->scE5x5() > 0.94) || (gsfiter->scE1x5()/gsfiter->scE5x5() > 0.83) ; break ;
    case kEndcap: result = true ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}

bool HEEPCut_41_isolEMHadDepth1::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  float gsf_ecaliso  = gsfiter->dr03EcalRecHitSumEt() ;
  float gsf_hcaliso1 = gsfiter->dr03HcalDepth1TowerSumEt() ;
  float gsf_gsfet   = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel:
      result = (gsf_ecaliso+gsf_hcaliso1) < (2.+0.03*gsf_gsfet + rho_*EcalHcal1EffAreaBarrel_) ;
      break ;
    case kEndcap:
      if(gsf_gsfet<50.0){ result = (gsf_ecaliso+gsf_hcaliso1) <  2.5                       + rho_*EcalHcal1EffAreaEndcaps_   ; }
      else              { result = (gsf_ecaliso+gsf_hcaliso1) < (2.5+0.03*(gsf_gsfet-50.0) + rho_*EcalHcal1EffAreaEndcaps_ ) ; }
      break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
void HEEPCut_41_isolEMHadDepth1::setRho(float rho){ rho_ = rho ; }
void HEEPCut_41_isolEMHadDepth1::setEcalHcal1EffAreaBarrel (float EcalHcal1EffAreaBarrel ){ EcalHcal1EffAreaBarrel_  = EcalHcal1EffAreaBarrel  ; }
void HEEPCut_41_isolEMHadDepth1::setEcalHcal1EffAreaEndcaps(float EcalHcal1EffAreaEndcaps){ EcalHcal1EffAreaEndcaps_ = EcalHcal1EffAreaEndcaps ; }

bool HEEPCut_41_HcalIso2       ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = true ; break ;
    case kEndcap: result = true ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_IsolPtTrks     ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  float gsf_trackIso = gsfiter->dr03TkSumPt() ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (gsf_trackIso < 5.0) ; break ;
    case kEndcap: result = (gsf_trackIso < 5.0) ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_EcalDriven     ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  bool result = gsfiter->ecalDrivenSeed() ;
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_invalid        ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  bool result = true ;
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_missingHits    ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  bool result = (gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() <= 1) ;
  setStatus(result) ;
  return getStatus() ;
}
bool HEEPCut_41_conversion     ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  bool result = (gsfiter->convFlags() != 3) ;
  setStatus(result) ;
  return getStatus() ;
}

bool HEEPCut_41_dxyFirstPV     ::applyCut(reco::GsfElectronCollection::const_iterator gsfiter){
  float gsf_dxy_firstPVtx = gsfiter->gsfTrack()->dxy(firstPV_) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (fabs(gsf_dxy_firstPVtx) < 0.02) ; break ;
    case kEndcap: result = (fabs(gsf_dxy_firstPVtx) < 0.05) ; break ;
    default : break ;
  }
  setStatus(result) ;
  return getStatus() ;
}
void HEEPCut_41_dxyFirstPV     ::setFirstPV(math::XYZPoint PV){ firstPV_ = math::XYZPoint(PV) ; }

HEEPCutCollection::HEEPCutCollection(std::string name){
  name_ = name ;
  int lower = -0.5 ;
  int upper = 10.5 ;
  int nBins = 1+upper-lower ;
  TH1F* hBaseNElectrons = new TH1F("hBaseHEEPNElectrons", "", nBins, lower, upper) ;
  hBaseNElectrons->GetXaxis()->SetTitle("N(e)") ;
  hBaseNElectrons->GetYaxis()->SetTitle("events") ;
  hNElectrons_     = (TH1F*) hBaseNElectrons->Clone("hHEEPNElectrons"    ) ;
  hNElectronsPass_ = (TH1F*) hBaseNElectrons->Clone("hHEEPNElectronsPass") ;
}
HEEPCutCollection::~HEEPCutCollection(){};

void HEEPCutCollection::addCut(HEEPCutBase* cut){
  listOfCuts_.push_back(cut) ;
  cutTypes_.push_back(kCut) ;
}
void HEEPCutCollection::addCutCollection(HEEPCutCollection* collection){
  listOfCutCollections_.push_back(collection) ;
  cutTypes_.push_back(kCollection) ;
}
bool HEEPCutCollection::applyCuts(reco::GsfElectronCollection::const_iterator gsfiter, IIHEModuleHEEP* mod){
  cutIndex_        = 0 ;
  collectionIndex_ = 0 ;
  status_ = true ;
  for(unsigned int i=0 ; i<cutTypes_.size() ; ++i){
    switch(cutTypes_.at(i)){
      case kCut:{
        HEEPCutBase* cut = (HEEPCutBase*) listOfCuts_.at(cutIndex_) ;
        cut->applyCut(gsfiter) ;
        cut->store(mod) ;
        ++cutIndex_ ;
        status_ = (status_ && cut->getStatus()) ;
        break ;
      }
      case kCollection:{
        HEEPCutCollection* collection = listOfCutCollections_.at(collectionIndex_) ;
        collection->applyCuts(gsfiter, mod) ;
        ++collectionIndex_ ;
        status_ = (status_ && collection->getStatus()) ;
        break ;
      }
    }
  }
  return status_ ;
}
void HEEPCutCollection::config(IIHEModuleHEEP* mod){
  cutIndex_        = 0 ;
  collectionIndex_ = 0 ;
  for(unsigned int i=0 ; i<cutTypes_.size() ; ++i){
    switch(cutTypes_.at(i)){
      case kCut:{
        listOfCuts_.at(cutIndex_)->addBranch(mod) ;
        ++cutIndex_ ;
        break ;
      }
      case kCollection:{
        listOfCutCollections_.at(collectionIndex_)->config(mod) ;
        ++collectionIndex_ ;
        break ;
      }
    }
  }
}


