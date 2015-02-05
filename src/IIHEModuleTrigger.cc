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
  hlTriggerResultsTag_ = iConfig.getParameter<edm::InputTag>("TriggerResults") ;
  nEvents_ = 0 ;
  nWasRun_ = 0 ;
  nAccept_ = 0 ;
  nErrors_ = 0 ;
  
  std::string   photonTriggerMatchingsIn = iConfig.getUntrackedParameter<std::string>("triggerPhotonMatchings"   , "") ;
  std::string electronTriggerMatchingsIn = iConfig.getUntrackedParameter<std::string>("triggerElectronMatchings" , "") ;
  std::string     muonTriggerMatchingsIn = iConfig.getUntrackedParameter<std::string>("triggerMuonMatchings"     , "") ;
  std::string      tauTriggerMatchingsIn = iConfig.getUntrackedParameter<std::string>("triggerTauMatchings"      , "") ;
  std::string      jetTriggerMatchingsIn = iConfig.getUntrackedParameter<std::string>("triggerJetMatchings"      , "") ;
  
  std::vector<std::string>   photonTriggerMatchings = splitString(  photonTriggerMatchingsIn, ",") ;
  std::vector<std::string> electronTriggerMatchings = splitString(electronTriggerMatchingsIn, ",") ;
  std::vector<std::string>     muonTriggerMatchings = splitString(    muonTriggerMatchingsIn, ",") ;
  std::vector<std::string>      tauTriggerMatchings = splitString(     tauTriggerMatchingsIn, ",") ;
  std::vector<std::string>      jetTriggerMatchings = splitString(     jetTriggerMatchingsIn, ",") ;
  
  for(unsigned i=0 ; i<photonTriggerMatchings.size() ; ++i){
    triggerMatchings_.push_back(new TriggerMatchParameters(kHighLevel, kPhoton  , photonTriggerMatchings.at(i), "")) ;
  }
  for(unsigned i=0 ; i<electronTriggerMatchings.size() ; ++i){
    triggerMatchings_.push_back(new TriggerMatchParameters(kHighLevel, kElectron, electronTriggerMatchings.at(i), "")) ;
  }
  for(unsigned i=0 ; i<muonTriggerMatchings.size() ; ++i){
    triggerMatchings_.push_back(new TriggerMatchParameters(kHighLevel, kMuon    ,     muonTriggerMatchings.at(i), "")) ;
  }
  for(unsigned i=0 ; i<tauTriggerMatchings.size() ; ++i){
    triggerMatchings_.push_back(new TriggerMatchParameters(kHighLevel, kTau     ,      tauTriggerMatchings.at(i), "")) ;
  }
  for(unsigned i=0 ; i<jetTriggerMatchings.size() ; ++i){
    triggerMatchings_.push_back(new TriggerMatchParameters(kHighLevel, kJet     ,      jetTriggerMatchings.at(i), "")) ;
  }
  
  std::string triggersIn = iConfig.getUntrackedParameter<std::string>("triggers" , "") ;
  triggerNamesFromPSet_ = splitString(triggersIn, ",") ;
  
  includeSingleElectronTriggers_ = (triggersIn.find("singleElectron")!=std::string::npos) ;
  includeDoubleElectronTriggers_ = (triggersIn.find("doubleElectron")!=std::string::npos) ;
  includeTripleElectronTriggers_ = (triggersIn.find("tripleElectron")!=std::string::npos) ;
  includeSingleMuonTriggers_     = (triggersIn.find("singleMuon"    )!=std::string::npos) ;
  includeSingleMuonTriggers_     = (triggersIn.find("doubleMuon"    )!=std::string::npos) ;
  includeSingleElectronSingleMuonTriggers_ = (triggersIn.find("singleElectronSingleMuon")!=std::string::npos) ;
  includeSingleElectronDoubleMuonTriggers_ = (triggersIn.find("singleElectronDoubleMuon")!=std::string::npos) ;
  includeDoubleElectronSingleMuonTriggers_ = (triggersIn.find("doubleElectronSingleMuon")!=std::string::npos) ;
}
IIHEModuleTrigger::~IIHEModuleTrigger(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleTrigger::beginJob(){
}

bool IIHEModuleTrigger::addHLTrigger(HLTrigger* hlt){
  for(unsigned int i=0 ; i<HLTriggers_.size() ; ++i){
    if(HLTriggers_.at(i)->name()==hlt->name()){
      return false ;
    }
  }
  HLTriggers_.push_back(hlt) ;
  return true ;
}

bool IIHEModuleTrigger::addHLTrigger(HLTrigger* hlt, std::vector<std::string> filterNames, std::vector<TriggerMatchParameters*> TMPs){
  for(unsigned int i=0 ; i<HLTriggers_.size() ; ++i){
    if(HLTriggers_.at(i)->name()==hlt->name()){
      return false ;
    }
  }
  
  for(unsigned i=0 ; i<filterNames.size() ; ++i){
    for(unsigned j=0 ; j<TMPs.size() ; ++j){
      TriggerMatchParameters* TMPtmp = new TriggerMatchParameters(TMPs.at(j), filterNames.at(i)) ;
      hlt->addMatching(TMPtmp) ;
    }
  }
  
  HLTriggers_.push_back(hlt) ;
  return true ;
}

int IIHEModuleTrigger::addBranches(){
  int result = 0 ;
  IIHEAnalysis* analysis = parent_ ;
  for(unsigned int i=0 ; i<HLTriggers_.size() ; i++){
    result += HLTriggers_.at(i)->createBranches(analysis, nEvents_) ;
  }
  return result ;
}

// ------------ method called to for each event  ------------
void IIHEModuleTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Trigger information
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT") ;
  edm::Handle<trigger::TriggerEvent> trigEvent ; 
  iEvent.getByLabel(trigEventTag,trigEvent) ;
  
  // get hold of TriggerResults
  edm::Handle<TriggerResults> HLTR;
  iEvent.getByLabel(hlTriggerResultsTag_, HLTR) ;
  
  // Now fill the values
  IIHEAnalysis* analysis = parent_ ;
  for(unsigned int i=0 ; i<HLTriggers_.size() ; i++){
    HLTrigger* hlt = HLTriggers_.at(i) ;
    hlt->status(iEvent, iSetup, hltConfig_, HLTR, trigEvent, trigEventTag) ;
    hlt->fill(analysis, trigEvent, trigEventTag) ;
  }
  
  nEvents_++ ;
}

void IIHEModuleTrigger::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){
  bool changed = true ;
  
  if(hltConfig_.init(iRun, iSetup, hlTriggerResultsTag_.process(), changed)){
    if(changed){
      if(false) hltConfig_.dump("Modules") ;
      
      // Get the updated list of trigger names
      HLTNamesFromConfig_ = hltConfig_.triggerNames() ;
      for(unsigned int i=0 ; i<HLTNamesFromConfig_.size() ; ++i){
        std::string name = HLTNamesFromConfig_.at(i) ;
        // Attempt to add the trigger
        bool addThisTrigger = false ;
        
        HLTrigger* hlt = new HLTrigger(name) ;
        
        // First check to see if it's in the list of requested triggers
        if(hlt->isOnlySingleElectron()           && includeSingleElectronTriggers_          ) addThisTrigger = true ;
        if(hlt->isOnlyDoubleElectron()           && includeDoubleElectronTriggers_          ) addThisTrigger = true ;
        if(hlt->isOnlyTripleElectron()           && includeTripleElectronTriggers_          ) addThisTrigger = true ;
        if(hlt->isOnlySingleMuon()               && includeSingleMuonTriggers_              ) addThisTrigger = true ;
        if(hlt->isOnlyDoubleMuon()               && includeSingleMuonTriggers_              ) addThisTrigger = true ;
        if(hlt->isOnlyTripleMuon()               && includeSingleMuonTriggers_              ) addThisTrigger = true ;
        if(hlt->isOnlySingleElectronSingleMuon() && includeSingleElectronSingleMuonTriggers_) addThisTrigger = true ;
        if(hlt->isOnlySingleElectronDoubleMuon() && includeSingleElectronDoubleMuonTriggers_) addThisTrigger = true ;
        if(hlt->isOnlyDoubleElectronSingleMuon() && includeDoubleElectronSingleMuonTriggers_) addThisTrigger = true ;
        
        // Only loop over trigger names if we have to
        if(addThisTrigger==false){
          for(unsigned int j=0 ; j<triggerNamesFromPSet_.size() ; ++j){
            if(triggerNamesFromPSet_.at(j)==name){
              addThisTrigger = true ;
              break ;
            }
          }
        }
        
        if(addThisTrigger==false) continue ;
        
        std::vector<TriggerMatchParameters*> matchParameters ;
        for(unsigned int j=0 ; j<triggerMatchings_.size() ; ++j){
          if(name==triggerMatchings_.at(j)->triggerName()){
            matchParameters.push_back(triggerMatchings_.at(j)) ;
          }
        }
        if(matchParameters.size()>0){
          std::vector<std::string> moduleNames = hltConfig_.moduleLabels(i) ;
          std::vector<std::string> moduleNamesWithTags ;
          for(unsigned int j=0 ; j<moduleNames.size() ; ++j){
            if(hltConfig_.saveTags(moduleNames.at(j))){
              moduleNamesWithTags.push_back(moduleNames.at(j)) ;
            }
          }
          addHLTrigger(hlt, moduleNamesWithTags, matchParameters) ;
        }
        else{
          addHLTrigger(hlt) ;
        }
      }
      
      // Now we need to re-map the indices to the names, given that some new triggers may have been inserted to the menu
      for(unsigned int i=0 ; i<HLTriggers_.size() ; ++i){
        HLTriggers_.at(i)->beginRun(HLTNamesFromConfig_) ;
      }
      
      // Attempt to add branches
      addBranches() ;
      parent_->configureBranches() ;
    
      // Now reset things to 0
      nEvents_ = 0 ;
      nWasRun_ = 0 ;
      nAccept_ = 0 ;
      nErrors_ = 0 ;
    }
  }
  else{
    std::cout << "Failed to init hltConfig" << std::endl ;
  }
}

void IIHEModuleTrigger::beginEvent(){}
void IIHEModuleTrigger::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTrigger::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleTrigger);
