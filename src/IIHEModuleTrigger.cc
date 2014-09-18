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
  for(unsigned int i=0 ; i<L1TriggerNamesElectronIn.size() ; ++i){
    addL1TriggerElecton(L1TriggerNamesElectronIn.at(i)) ;
  }
  for(unsigned int i=0 ; i<HLTriggerNamesElectronIn.size() ; ++i){
    addHLTriggerElecton(HLTriggerNamesElectronIn.at(i).first, HLTriggerNamesElectronIn.at(i).second) ;
  }
  addBranches() ;
}

bool IIHEModuleTrigger::addL1TriggerElecton(std::string name){
  for(unsigned int i=0 ; i<L1TriggersElectron_.size() ; i++){
    if(name==L1TriggersElectron_.at(i)->name()){
      return false ;
    }
  }
  L1TriggersElectron_.push_back(new L1Trigger(name, branchPrefixElectronMatch_)) ;
  return true ;
}

bool IIHEModuleTrigger::addHLTriggerElecton(std::string name, float DeltaRCut){
  for(unsigned int i=0 ; i<HLTriggersElectron_ .size(); i++){
    if(name==HLTriggersElectron_.at(i)->name()){
      return false ;
    }
  }
  HLTriggersElectron_.push_back(new HLTrigger(name, branchPrefixElectronMatch_, DeltaRCut)) ;
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
}

// ------------ method called to for each event  ------------
void IIHEModuleTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  
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
    for(unsigned int i=0 ; i<L1TriggersElectron_.size() ; i++){
      store(L1TriggersElectron_.at(i)->branchName(), L1TriggersElectron_.at(i)->matchElectron(trigEvent, gsfiter)) ;
    }
    for(unsigned int i=0 ; i<HLTriggersElectron_.size() ; i++){
      store(HLTriggersElectron_.at(i)->branchName(), HLTriggersElectron_.at(i)->matchElectron(trigEvent, gsfiter)) ;
    }
  }
  
}

void IIHEModuleTrigger::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}

void IIHEModuleTrigger::beginEvent(){}
void IIHEModuleTrigger::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTrigger::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleTrigger);
