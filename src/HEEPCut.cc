#include "UserCode/IIHETree/interface/HEEPCut.h"

#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

//////////////////////////////////////////////////////////////////////////////////////////
//                                      Base class                                      //
//////////////////////////////////////////////////////////////////////////////////////////
HEEPCutBase::HEEPCutBase(std::string name, IIHEModuleHEEP* mod){
  name_ = name ;
  index_ = -1 ;
  nPass_ = 0 ;
  nPassCumulative_ = 0 ;
  value_ = -999 ;
  reset() ;
  parent_mod_ = mod ;
  branchName_n_           = name_ + "_n" ;
  branchName_nCumulative_ = name_ + "_nCumulative" ;
  branchName_value_       = name_ + "_value" ;
}
HEEPCutBase::~HEEPCutBase(){}
void HEEPCutBase::setStatus(bool value, bool cumulativeSuccess){
  status_ = value ;
  if(status_) nPass_++ ;
  if(status_ && cumulativeSuccess) nPassCumulative_++ ;
}
bool HEEPCutBase::getStatus(){ return status_ ; }
void HEEPCutBase::reset(){ status_ = true ; }
bool HEEPCutBase::addBranches(){
  bool success = true ;
  success = (success && parent_mod_->addBranch(name_                  , kVectorBool )) ;
  success = (success && parent_mod_->addBranch(branchName_n_          , kInt        )) ;
  success = (success && parent_mod_->addBranch(branchName_nCumulative_, kInt        )) ;
  success = (success && parent_mod_->addBranch(branchName_value_      , kVectorFloat)) ;
  return success ;
}
void HEEPCutBase::store(){
  parent_mod_->store(name_, status_) ;
  parent_mod_->store(branchName_value_, value_) ;
}
bool HEEPCutBase::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){ return status_ ; }
void HEEPCutBase::beginEvent(){
  nPass_           = 0 ;
  nPassCumulative_ = 0 ;
}
void HEEPCutBase::endEvent(){
  parent_mod_->store(branchName_n_          , nPass_          ) ;
  parent_mod_->store(branchName_nCumulative_, nPassCumulative_) ;
}

bool HEEPCutBase::isBarrel(reco::GsfElectron* gsfiter){
  return (fabs(gsfiter->superCluster()->eta()) < barrelEtaUpper_ ) ;
}
bool HEEPCutBase::isEndcap(reco::GsfElectron* gsfiter){
  float eta = gsfiter->superCluster()->eta() ;
  return (fabs(eta) > endcapEtaLower_ && fabs(eta) < endcapEtaUpper_) ;
}
int  HEEPCutBase::detectorRegion(reco::GsfElectron* gsfiter){
  if(isBarrel(gsfiter)) return kBarrel ;
  if(isEndcap(gsfiter)) return kEndcap ;
  return kNone ;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                    Base HEEP cuts                                    //
//////////////////////////////////////////////////////////////////////////////////////////
HEEPCut_Et::HEEPCut_Et(std::string name, IIHEModuleHEEP* mod, float thresholdBarrel, float thresholdEndcap):
  HEEPCutBase(name, mod){
  thresholdBarrel_ = thresholdBarrel ;
  thresholdEndcap_ = thresholdEndcap ;
} ;
bool HEEPCut_Et::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->caloEnergy()*sin(gsfiter->p4().theta())) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (value() > thresholdBarrel_) ; break ;
    case kEndcap: result = (value() > thresholdEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

HEEPCut_eta::HEEPCut_eta(std::string name, IIHEModuleHEEP* mod): HEEPCutBase(name, mod){} ;
bool HEEPCut_eta::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->superCluster()->eta()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = (region==kBarrel || region==kEndcap) ;
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

HEEPCut_EcalDriven::HEEPCut_EcalDriven(std::string name, IIHEModuleHEEP* mod): HEEPCutBase(name, mod){} ;
bool HEEPCut_EcalDriven::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  bool result = gsfiter->ecalDrivenSeed() ;
  setStatus(result, cumulativeSuccess) ;
  setValue((float) result) ;
  return getStatus() ;
}

HEEPCut_dPhiIn::HEEPCut_dPhiIn(std::string name, IIHEModuleHEEP* mod, float thresholdBarrel, float thresholdEndcap): 
  HEEPCutBase(name, mod){
  thresholdBarrel_ = thresholdBarrel ;
  thresholdEndcap_ = thresholdEndcap ;
} ;
bool HEEPCut_dPhiIn::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->deltaPhiSuperClusterTrackAtVtx()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (fabs(value()) < thresholdBarrel_) ; break ;
    case kEndcap: result = (fabs(value()) < thresholdEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

HEEPCut_SigmaIetaIeta::HEEPCut_SigmaIetaIeta(std::string name, IIHEModuleHEEP* mod, float threshold):
  HEEPCutBase(name, mod){
  threshold_ = threshold ;
} ;
bool HEEPCut_SigmaIetaIeta::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->sigmaIetaIeta()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = true ; break ;
    case kEndcap: result = (value() < threshold_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

// These cuts go together
HEEPCut_E1x5OverE5x5::HEEPCut_E1x5OverE5x5(std::string name, IIHEModuleHEEP* mod): HEEPCutBase(name, mod){} ;
HEEPCut_E2x5OverE5x5::HEEPCut_E2x5OverE5x5(std::string name, IIHEModuleHEEP* mod, float thresholdE1x5, float thresholdE2x5):
  HEEPCutBase(name, mod){
  thresholdE1x5_ = thresholdE1x5 ;
  thresholdE2x5_ = thresholdE2x5 ;
} ;
bool HEEPCut_E1x5OverE5x5::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->scE1x5()/gsfiter->scE5x5()) ;
  bool result = true ;
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}
bool HEEPCut_E2x5OverE5x5::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->scE2x5Max()/gsfiter->scE5x5()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (gsfiter->scE2x5Max()/gsfiter->scE5x5() > thresholdE2x5_) || (gsfiter->scE1x5()/gsfiter->scE5x5() > thresholdE1x5_) ; break ;
    case kEndcap: result = true ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}


HEEPCut_isolEMHadDepth1::HEEPCut_isolEMHadDepth1(std::string name, IIHEModuleHEEP* mod, float constantTermBarrel, float constantTermEndcapLowEt, float constantTermEndcapHighEt, float linearTermBarrel, float linearTermEndcap, float offsetTermEndcap):
  HEEPCutBase(name, mod){
  constantTermBarrel_       = constantTermBarrel       ;
  constantTermEndcapLowEt_  = constantTermEndcapLowEt  ;
  constantTermEndcapHighEt_ = constantTermEndcapHighEt ;
  linearTermBarrel_         = linearTermBarrel         ;
  linearTermEndcap_         = linearTermEndcap         ;
  offsetTermEndcap_         = offsetTermEndcap         ;
} ;
bool HEEPCut_isolEMHadDepth1::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  float gsf_ecaliso  = gsfiter->dr03EcalRecHitSumEt() ;
  float gsf_hcaliso1 = gsfiter->dr03HcalDepth1TowerSumEt() ;
  float gsf_gsfet   = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
  setValue(gsf_ecaliso+gsf_hcaliso1) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel:
      result = (value()) < (constantTermBarrel_ + linearTermBarrel_*gsf_gsfet + rho_*EcalHcal1EffAreaBarrel_) ;
      break ;
    case kEndcap:
      if(gsf_gsfet<offsetTermEndcap_){ result = value() <  constantTermEndcapLowEt_                                                    + rho_*EcalHcal1EffAreaEndcaps_   ; }
      else                           { result = value() < (constantTermEndcapHighEt_ + linearTermEndcap_*(gsf_gsfet-offsetTermEndcap_) + rho_*EcalHcal1EffAreaEndcaps_ ) ; }
      break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}
void HEEPCut_isolEMHadDepth1::setRho(float rho){ rho_ = rho ; }
void HEEPCut_isolEMHadDepth1::setEcalHcal1EffAreaBarrel (float EcalHcal1EffAreaBarrel ){ EcalHcal1EffAreaBarrel_  = EcalHcal1EffAreaBarrel  ; }
void HEEPCut_isolEMHadDepth1::setEcalHcal1EffAreaEndcaps(float EcalHcal1EffAreaEndcaps){ EcalHcal1EffAreaEndcaps_ = EcalHcal1EffAreaEndcaps ; }

HEEPCut_IsolPtTrks::HEEPCut_IsolPtTrks(std::string name, IIHEModuleHEEP* mod, float thresholdBarrel, float thresholdEndcap):
  HEEPCutBase(name, mod){
  thresholdBarrel_ = thresholdBarrel ;
  thresholdEndcap_ = thresholdEndcap ;
} ;
bool HEEPCut_IsolPtTrks::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->dr03TkSumPt()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (value() < thresholdBarrel_) ; break ;
    case kEndcap: result = (value() < thresholdEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

HEEPCut_missingHits::HEEPCut_missingHits(std::string name, IIHEModuleHEEP* mod, float threshold):
  HEEPCutBase(name, mod){
  threshold_ = threshold ;
} ;
bool HEEPCut_missingHits::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue((float) (gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits())) ;
  bool result = (value() <= threshold_+0.5) ;
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

HEEPCut_dxyFirstPV ::HEEPCut_dxyFirstPV(std::string name, IIHEModuleHEEP* mod, float thresholdBarrel, float thresholdEndcap):
  HEEPCutBase(name, mod){
  thresholdBarrel_ = thresholdBarrel ;
  thresholdEndcap_ = thresholdEndcap ;
} ;
bool HEEPCut_dxyFirstPV::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->gsfTrack()->dxy(firstPV_)) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (fabs(value()) < thresholdBarrel_) ; break ;
    case kEndcap: result = (fabs(value()) < thresholdEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}
void HEEPCut_dxyFirstPV::setFirstPV(math::XYZPoint* PV){ firstPV_ = math::XYZPoint(*PV) ; }


//////////////////////////////////////////////////////////////////////////////////////////
//                                     HEEP 4.1 cuts                                    //
//////////////////////////////////////////////////////////////////////////////////////////

HEEPCut_41_dEtaIn::HEEPCut_41_dEtaIn(std::string name, IIHEModuleHEEP* mod, float thresholdBarrel, float thresholdEndcap):
  HEEPCutBase(name, mod){
  thresholdBarrel_ = thresholdBarrel ;
  thresholdEndcap_ = thresholdEndcap ;
} ;

HEEPCut_41_HOverE::HEEPCut_41_HOverE(std::string name, IIHEModuleHEEP* mod, float thresholdBarrel, float thresholdEndcap): 
  HEEPCutBase(name, mod){
  thresholdBarrel_ = thresholdBarrel ;
  thresholdEndcap_ = thresholdEndcap ;
} ;

bool HEEPCut_41_dEtaIn::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->deltaEtaSuperClusterTrackAtVtx()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (fabs(value()) < thresholdBarrel_) ; break ;
    case kEndcap: result = (fabs(value()) < thresholdEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

bool HEEPCut_41_HOverE::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->hadronicOverEm()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (value() < thresholdBarrel_) ; break ;
    case kEndcap: result = (value() < thresholdEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                 HEEP 5.0  startup cuts                               //
//////////////////////////////////////////////////////////////////////////////////////////
HEEPCut_50_50ns_dEtaIn::HEEPCut_50_50ns_dEtaIn(std::string name, IIHEModuleHEEP* mod, float constantTermBarrel, float linearTermBarrel, float cutoffTermBarrel, float thresholdEndcap):
  HEEPCutBase(name, mod){
  constantTermBarrel_ = constantTermBarrel ;
  linearTermBarrel_   = linearTermBarrel   ;
  cutoffTermBarrel_   = cutoffTermBarrel   ;
  thresholdEndcap_    = thresholdEndcap    ;
} ;
HEEPCut_50_50ns_HOverE::HEEPCut_50_50ns_HOverE(std::string name, IIHEModuleHEEP* mod, float reciprocalTermBarrel, float reciprocalTermEndcap, float constantTermBarrel, float constantTermEndcap):
  HEEPCutBase(name, mod){
  reciprocalTermBarrel_ = reciprocalTermBarrel_ ;
  reciprocalTermEndcap_ = reciprocalTermEndcap_ ;
  constantTermBarrel_   = constantTermBarrel_   ;
  constantTermEndcap_   = constantTermEndcap_   ;
} ;

bool HEEPCut_50_50ns_dEtaIn::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->deltaEtaSuperClusterTrackAtVtx()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel:{
      float Et  = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
      float threshold = std::max(constantTermBarrel_ - linearTermBarrel_*Et, cutoffTermBarrel_) ;
      result = (fabs(value()) < threshold) ;
      break ;
    }
    case kEndcap:{
      result = (fabs(value()) < thresholdEndcap_) ;
      break ;
    }
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}
bool HEEPCut_50_50ns_HOverE::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->hadronicOverEm()) ;
  float E = gsfiter->caloEnergy() ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (value() < reciprocalTermBarrel_/E + constantTermBarrel_) ; break ;
    case kEndcap: result = (value() < reciprocalTermEndcap_/E + constantTermEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                   HEEP 5.0 25ns cuts                                 //
//////////////////////////////////////////////////////////////////////////////////////////
HEEPCut_50_25ns_dEtaIn::HEEPCut_50_25ns_dEtaIn(std::string name, IIHEModuleHEEP* mod, float constantTermBarrel, float linearTermBarrel, float cutoffTermBarrel, float constantTermEndcap, float linearTermEndcap, float cutoffTermEndcap):
  HEEPCutBase(name, mod){
  constantTermBarrel_ = constantTermBarrel ;
  linearTermBarrel_   = linearTermBarrel   ;
  cutoffTermBarrel_   = cutoffTermBarrel   ;
  constantTermEndcap_ = constantTermEndcap ;
  linearTermEndcap_   = linearTermEndcap   ;
  cutoffTermEndcap_   = cutoffTermEndcap   ;
} ;
HEEPCut_50_25ns_HOverE::HEEPCut_50_25ns_HOverE(std::string name, IIHEModuleHEEP* mod, float reciprocalTermBarrel, float constantTermBarrel, float reciprocalTermEndcap, float constantTermEndcap):
  HEEPCutBase(name, mod){
  reciprocalTermBarrel_ = reciprocalTermBarrel_ ;
  constantTermBarrel_   = constantTermBarrel_   ;
  reciprocalTermEndcap_ = reciprocalTermEndcap_ ;
  constantTermEndcap_   = constantTermEndcap_   ;
} ;

bool HEEPCut_50_25ns_dEtaIn::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->deltaEtaSuperClusterTrackAtVtx()) ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  float Et = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
  switch(region){
    case kBarrel:{
      float threshold = std::max(constantTermBarrel_ - constantTermBarrel_*Et, cutoffTermBarrel_) ;
      result = (fabs(value()) < threshold) ;
      break ;
    }
    case kEndcap:{
      float threshold = std::max(constantTermBarrel_ - constantTermBarrel_*Et, cutoffTermBarrel_) ;
      result = (fabs(value()) < threshold) ;
      break ;
    }
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}
bool HEEPCut_50_25ns_HOverE::applyCut(reco::GsfElectron* gsfiter, bool cumulativeSuccess){
  setValue(gsfiter->hadronicOverEm()) ;
  float E = gsfiter->caloEnergy() ;
  int region = detectorRegion(gsfiter) ;
  bool result = true ;
  switch(region){
    case kBarrel: result = (value() < reciprocalTermBarrel_/E + constantTermBarrel_) ; break ;
    case kEndcap: result = (value() < reciprocalTermEndcap_/E + constantTermEndcap_) ; break ;
    default : break ;
  }
  setStatus(result, cumulativeSuccess) ;
  return getStatus() ;
}

//////////////////////////////////////////////////////////////////////////////////////////
//                                   Cut collections                                    //
//////////////////////////////////////////////////////////////////////////////////////////
HEEPCutCollection::HEEPCutCollection(std::string name, IIHEModuleHEEP* mod){
  name_         = name ;
  branchName_n_ = name + "_n" ;
  parent_mod_   = mod ;
  nPass_        = -1 ;
}
HEEPCutCollection::~HEEPCutCollection(){};
void HEEPCutCollection::config(float barrelEtaUpper, float endcapEtaLower, float endcapEtaUpper){
  cutIndex_        = 0 ;
  collectionIndex_ = 0 ;
  for(unsigned int i=0 ; i<cutTypes_.size() ; ++i){
    switch(cutTypes_.at(i)){
      case kCut:{
        HEEPCutBase* cut = listOfCuts_.at(cutIndex_) ;
        cut->addBranches() ;
        cut->setBarrelLimits(barrelEtaUpper) ;
        cut->setEndcapLimits(endcapEtaLower, endcapEtaUpper) ;
        ++cutIndex_ ;
        break ;
      }
      case kCollection:{
        listOfCutCollections_.at(collectionIndex_)->config(barrelEtaUpper, endcapEtaLower, endcapEtaUpper) ;
        ++collectionIndex_ ;
        break ;
      }
    }
  }
  parent_mod_->addBranch(name_, kVectorBool) ;
  parent_mod_->addBranch(branchName_n_, kInt) ;
}

void HEEPCutCollection::beginEvent(){
  nPass_           = 0 ;
  cutIndex_        = 0 ;
  collectionIndex_ = 0 ;
  for(unsigned int i=0 ; i<cutTypes_.size() ; ++i){
    switch(cutTypes_.at(i)){
      case kCut:{
        listOfCuts_.at(cutIndex_)->beginEvent() ;
        ++cutIndex_ ;
        break ;
      }
      case kCollection:{
        listOfCutCollections_.at(collectionIndex_)->beginEvent() ;
        ++collectionIndex_ ;
        break ;
      }
    }
  }
}
void HEEPCutCollection::endEvent()  {
  parent_mod_->store(branchName_n_, nPass_) ;
  cutIndex_        = 0 ;
  collectionIndex_ = 0 ;
  for(unsigned int i=0 ; i<cutTypes_.size() ; ++i){
    switch(cutTypes_.at(i)){
      case kCut:{
        listOfCuts_.at(cutIndex_)->endEvent() ;
        ++cutIndex_ ;
        break ;
      }
      case kCollection:{
        listOfCutCollections_.at(collectionIndex_)->endEvent() ;
        ++collectionIndex_ ;
        break ;
      }
    }
  }
}

void HEEPCutCollection::addCut(HEEPCutBase* cut){
  listOfCuts_.push_back(cut) ;
  cutTypes_.push_back(kCut) ;
}
void HEEPCutCollection::addCutCollection(HEEPCutCollection* collection){
  listOfCutCollections_.push_back(collection) ;
  cutTypes_.push_back(kCollection) ;
}
bool HEEPCutCollection::applyCuts(reco::GsfElectron* gsfiter){
  return applyCuts(gsfiter, true) ;
}
bool HEEPCutCollection::applyCuts(reco::GsfElectron* gsfiter, bool status_in){
  cutIndex_        = 0 ;
  collectionIndex_ = 0 ;
  status_ = status_in ;
  for(unsigned int i=0 ; i<cutTypes_.size() ; ++i){
    switch(cutTypes_.at(i)){
      case kCut:{
        HEEPCutBase* cut = (HEEPCutBase*) listOfCuts_.at(cutIndex_) ;
        cut->applyCut(gsfiter, status_) ;
        cut->store() ;
        ++cutIndex_ ;
        status_ = (status_ && cut->getStatus()) ;
        break ;
      }
      case kCollection:{
        HEEPCutCollection* collection = listOfCutCollections_.at(collectionIndex_) ;
        collection->applyCuts(gsfiter, status_) ;
        ++collectionIndex_ ;
        status_ = (status_ && collection->getStatus()) ;
        break ;
      }
    }
  }
  parent_mod_->store(name_, status_) ;
  if(status_) nPass_++ ;
  return status_ ;
}


