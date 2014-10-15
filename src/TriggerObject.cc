#include "UserCode/IIHETree/interface/TriggerObject.h"

L1Trigger::L1Trigger(std::string name, std::string prefix):
filterIndex_(-1)
,barrelEnd_(1.4791)
,regionEtaSizeEB_(0.522)
,regionEtaSizeEE_(1.0)
,regionPhiSize_(1.044){
  name_ = name ;
  branchName_ = prefix + name_ ;
}
L1Trigger::~L1Trigger(){}

int L1Trigger::setFilterIndex(edm::Handle<trigger::TriggerEvent> trigEvent, edm::InputTag trigEventTag){
  filterIndex_ = trigEvent->filterIndex(edm::InputTag(name_,"",trigEventTag.process())) ;
  return filterIndex_ ;
}

bool L1Trigger::matchElectron(edm::Handle<trigger::TriggerEvent> trigEvent, reco::GsfElectron* gsf){
  float eta = gsf->eta() ;
  float phi = gsf->phi() ;
  return matchObject(trigEvent, eta, phi) ;
}

bool L1Trigger::matchMuon(edm::Handle<trigger::TriggerEvent> trigEvent, reco::Muon* muon){
  float eta = muon->eta() ;
  float phi = muon->phi() ;
  return matchObject(trigEvent, eta, phi) ;
}

bool L1Trigger::matchObject(edm::Handle<trigger::TriggerEvent> trigEvent, float eta, float phi){
  // Careful that L1 triggers only have discrete eta phi. Need to be extremely loose. 
  // See here: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SHarper/SHNtupliser/src/SHTrigInfo.cc?revision=1.5&view=markup&pathrev=HEAD
  // It is important to specify the right HLT process for the filter, not doing this is a common bug
  if(filterIndex_<0) return false ;
  if(filterIndex_<trigEvent->sizeFilters()){ 
    const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex_) ;
    
    const trigger::TriggerObjectCollection& trigObjColl(trigEvent->getObjects()) ;
    
    // Now loop of the trigger objects passing filter
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin() ; keyIt!=trigKeys.end() ; ++keyIt){ 
      
      const trigger::TriggerObject& obj = trigObjColl[*keyIt] ;
      
      float objeta = obj.eta() ;
      float objphi = obj.phi() ;
      
      double etaBinLow  = 0.0 ;
      double etaBinHigh = 0.0 ;
      
      if(fabs(objeta)<barrelEnd_){
        etaBinLow  = objeta    - regionEtaSizeEB_/2.0 ;
        etaBinHigh = etaBinLow + regionEtaSizeEB_     ;
      }
      else{
        etaBinLow  = objeta    - regionEtaSizeEE_/2.0 ;
        etaBinHigh = etaBinLow + regionEtaSizeEE_     ;
      }
      
      if(eta<etaBinHigh && eta>etaBinLow){
        return true ;
      }
      float dPhi = reco::deltaPhi(phi,objphi) ;
      if(eta<etaBinHigh && eta>etaBinLow &&  dPhi<regionPhiSize_/2.0){
        return true ;
      }
    }
  }
  return false ;
}

HLTrigger::HLTrigger(std::string name, std::string prefix, float DeltaRCut=0.5):
filterIndex_(-1){
  name_ = name ;
  branchName_ = prefix + name_ ;
  DeltaRCut_ = DeltaRCut ;
}
HLTrigger::~HLTrigger(){}

int HLTrigger::setFilterIndex(edm::Handle<trigger::TriggerEvent> trigEvent, edm::InputTag trigEventTag){
  filterIndex_ = trigEvent->filterIndex(edm::InputTag(name_,"",trigEventTag.process())) ;
  return filterIndex_ ;
}

bool HLTrigger::matchElectron(edm::Handle<trigger::TriggerEvent> trigEvent, reco::GsfElectron* gsf){
  if(filterIndex_<0) return false ;
  float eta = gsf->eta() ;
  float phi = gsf->phi() ;
  return matchObject(trigEvent, eta, phi) ;
}

bool HLTrigger::matchMuon(edm::Handle<trigger::TriggerEvent> trigEvent, reco::Muon* muon){
  if(filterIndex_<0) return false ;
  float eta = muon->eta() ;
  float phi = muon->phi() ;
  return matchObject(trigEvent, eta, phi) ;
}

bool HLTrigger::matchObject(edm::Handle<trigger::TriggerEvent> trigEvent, float eta, float phi){
  if(filterIndex_<0) return false ;
  if(filterIndex_<trigEvent->sizeFilters()){ 
    const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex_) ; 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects()) ;
    
    // Now loop over the trigger objects passing filter
    for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); ++keyIt) { 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt] ;
      if(deltaR(eta,phi,obj.eta(),obj.phi())<DeltaRCut_){
        return true ;
      }
    }
  }
  return false ;
}

