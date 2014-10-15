#include "UserCode/IIHETree/interface/IIHEModuleTrigger.h"
#include "UserCode/IIHETree/interface/TriggerObject.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
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

IIHEModuleTrigger::IIHEModuleTrigger(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
  branchPrefixElectronMatch_ = "gsf_triggerMatch_" ;
}
IIHEModuleTrigger::~IIHEModuleTrigger(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleTrigger::beginJob(){
  std::vector<std::string>                   L1TriggerNamesElectronIn = parent_->getTriggerL1FilterNamesElectron () ;
  std::vector<std::pair<std::string,float> > HLTriggerNamesElectronIn = parent_->getTriggerHLTFilterNamesElectron() ;
  std::vector<std::string>                   L1TriggerNamesMuonIn     = parent_->getTriggerL1FilterNamesMuon     () ;
  std::vector<std::pair<std::string,float> > HLTriggerNamesMuonIn     = parent_->getTriggerHLTFilterNamesMuon    () ;
  for(unsigned int i=0 ; i<L1TriggerNamesElectronIn.size() ; ++i){
    addL1TriggerElectron(L1TriggerNamesElectronIn.at(i)) ;
  }
  for(unsigned int i=0 ; i<HLTriggerNamesElectronIn.size() ; ++i){
    addHLTriggerElectron(HLTriggerNamesElectronIn.at(i).first, HLTriggerNamesElectronIn.at(i).second) ;
  }
  for(unsigned int i=0 ; i<L1TriggerNamesMuonIn.size() ; ++i){
    addL1TriggerMuon(L1TriggerNamesMuonIn.at(i)) ;
  }
  for(unsigned int i=0 ; i<HLTriggerNamesMuonIn.size() ; ++i){
    addHLTriggerMuon(HLTriggerNamesMuonIn.at(i).first, HLTriggerNamesMuonIn.at(i).second) ;
  }
  addBranches() ;
}

bool IIHEModuleTrigger::addL1TriggerElectron(std::string name){
  for(unsigned int i=0 ; i<L1TriggersElectron_.size() ; i++){
    if(name==L1TriggersElectron_.at(i)->name()){
      return false ;
    }
  }
  L1TriggersElectron_.push_back(new L1Trigger(name, branchPrefixElectronMatch_)) ;
  return true ;
}
bool IIHEModuleTrigger::addHLTriggerElectron(std::string name, float DeltaRCut){
  for(unsigned int i=0 ; i<HLTriggersElectron_ .size(); i++){
    if(name==HLTriggersElectron_.at(i)->name()){
      return false ;
    }
  }
  HLTriggersElectron_.push_back(new HLTrigger(name, branchPrefixElectronMatch_, DeltaRCut)) ;
  return true ;
}

bool IIHEModuleTrigger::addL1TriggerMuon(std::string name){
  for(unsigned int i=0 ; i<L1TriggersMuon_.size() ; i++){
    if(name==L1TriggersMuon_.at(i)->name()){
      return false ;
    }
  }
  L1TriggersMuon_.push_back(new L1Trigger(name, branchPrefixMuonMatch_)) ;
  return true ;
}
bool IIHEModuleTrigger::addHLTriggerMuon(std::string name, float DeltaRCut){
  for(unsigned int i=0 ; i<HLTriggersMuon_ .size(); i++){
    if(name==HLTriggersMuon_.at(i)->name()){
      return false ;
    }
  }
  HLTriggersMuon_.push_back(new HLTrigger(name, branchPrefixMuonMatch_, DeltaRCut)) ;
  return true ;
}

void IIHEModuleTrigger::addBranches(){
  setBranchType(kVectorBool) ;
  for(unsigned int i=0 ; i<L1TriggersElectron_.size() ; i++){
    addBranch(L1TriggersElectron_.at(i)->branchName()) ;
  }
  for(unsigned int i=0 ; i<HLTriggersElectron_.size() ; i++){
    addBranch(HLTriggersElectron_.at(i)->branchName()) ;
  }
  for(unsigned int i=0 ; i<L1TriggersMuon_.size()     ; i++){
    addBranch(L1TriggersMuon_.at(i)->branchName()) ;
  }
  for(unsigned int i=0 ; i<HLTriggersMuon_.size()     ; i++){
    addBranch(HLTriggersMuon_.at(i)->branchName()) ;
  }
}

// ------------ method called to for each event  ------------
void IIHEModuleTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  reco::MuonCollection        muons     = parent_->getMuonCollection()     ;
  
  // Trigger information
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
  edm::Handle<trigger::TriggerEvent> trigEvent ; 
  iEvent.getByLabel(trigEventTag,trigEvent) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                  Trigger matching                                  //
  ////////////////////////////////////////////////////////////////////////////////////////
  // Electrons
  for(unsigned int i=0 ; i<L1TriggersElectron_.size() ; i++){
    L1TriggersElectron_.at(i)->setFilterIndex(trigEvent, trigEventTag) ;
  }
  for(unsigned int i=0 ; i<HLTriggersElectron_.size() ; i++){
    HLTriggersElectron_.at(i)->setFilterIndex(trigEvent, trigEventTag) ;
  }
  
  for(reco::GsfElectronCollection::const_iterator gsfiter=electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
    reco::GsfElectron* gsf = (reco::GsfElectron*) &*gsfiter ;
    for(unsigned int i=0 ; i<L1TriggersElectron_.size() ; i++){
      store(L1TriggersElectron_.at(i)->branchName(), L1TriggersElectron_.at(i)->matchElectron(trigEvent, gsf)) ;
    }
    for(unsigned int i=0 ; i<HLTriggersElectron_.size() ; i++){
      store(HLTriggersElectron_.at(i)->branchName(), HLTriggersElectron_.at(i)->matchElectron(trigEvent, gsf)) ;
    }
  }
  
  // Muons
  for(unsigned int i=0 ; i<L1TriggersMuon_.size() ; i++){
    L1TriggersMuon_.at(i)->setFilterIndex(trigEvent, trigEventTag) ;
  }
  for(unsigned int i=0 ; i<HLTriggersMuon_.size() ; i++){
    HLTriggersMuon_.at(i)->setFilterIndex(trigEvent, trigEventTag) ;
  }
  
  for(reco::MuonCollection::const_iterator muiter=muons.begin() ; muiter!=muons.end() ; ++muiter){
    reco::Muon* muon = (reco::Muon*) &*muiter ;
    for(unsigned int i=0 ; i<L1TriggersMuon_.size() ; i++){
      store(L1TriggersMuon_.at(i)->branchName(), L1TriggersMuon_.at(i)->matchMuon(trigEvent, muon)) ;
    }
    for(unsigned int i=0 ; i<HLTriggersMuon_.size() ; i++){
      store(HLTriggersMuon_.at(i)->branchName(), HLTriggersMuon_.at(i)->matchMuon(trigEvent, muon)) ;
    }
  }
}

void IIHEModuleTrigger::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}

void IIHEModuleTrigger::beginEvent(){}
void IIHEModuleTrigger::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTrigger::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleTrigger);
