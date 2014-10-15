// System includes
#include <iostream>
#include <TMath.h>
#include <vector>

// CMSSW includes
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Common/interface/TriggerNames.h"

// IIHE includes
#include "UserCode/IIHETree/interface/IIHEAnalysis.h"

#include "UserCode/IIHETree/interface/BranchWrapper.h"
#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "UserCode/IIHETree/interface/IIHEModuleEvent.h"
#include "UserCode/IIHETree/interface/IIHEModuleVertex.h"
#include "UserCode/IIHETree/interface/IIHEModuleSuperCluster.h"
#include "UserCode/IIHETree/interface/IIHEModulePhoton.h"
#include "UserCode/IIHETree/interface/IIHEModuleGedGsfElectron.h"
#include "UserCode/IIHETree/interface/IIHEModuleMuon.h"
#include "UserCode/IIHETree/interface/IIHEModuleHEEP.h"
#include "UserCode/IIHETree/interface/IIHEModuleMET.h"
#include "UserCode/IIHETree/interface/IIHEModuleMCTruth.h"
#include "UserCode/IIHETree/interface/IIHEModuleTrigger.h"

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEAnalysis::IIHEAnalysis(const edm::ParameterSet& iConfig){
  currentVarType_ = -1 ;
  debug_     = iConfig.getParameter<bool  >("debug"    ) ;
  git_hash_  = iConfig.getParameter<string>("git_hash" ) ;
  globalTag_ = iConfig.getParameter<string>("globalTag") ;
  beamSpotLabel_      = consumes<BeamSpot>(iConfig.getParameter<InputTag>("beamSpot")) ;
  primaryVertexLabel_ = iConfig.getParameter<edm::InputTag>("primaryVertex") ;
  
  photonCollectionLabel_   = iConfig.getParameter<edm::InputTag>("photonCollection"       ) ;
  electronCollectionLabel_ = iConfig.getParameter<edm::InputTag>("electronCollection"     ) ;
  muonCollectionLabel_     = iConfig.getParameter<edm::InputTag>("muonCollection"         ) ;
  
  IIHEModuleTrigger* mod_trigger = new IIHEModuleTrigger(iConfig) ;
  
  childModules_.push_back(new IIHEModuleEvent(iConfig)         ) ;
  childModules_.push_back(new IIHEModuleVertex(iConfig)        ) ;
  childModules_.push_back(new IIHEModuleSuperCluster(iConfig)  ) ;
  childModules_.push_back(new IIHEModulePhoton(iConfig)        ) ;
  childModules_.push_back(new IIHEModuleGedGsfElectron(iConfig)) ;
  childModules_.push_back(new IIHEModuleMuon(iConfig)          ) ;
  childModules_.push_back(new IIHEModuleMET(iConfig)           ) ;
  childModules_.push_back(new IIHEModuleHEEP(iConfig)          ) ;
  childModules_.push_back(new IIHEModuleMCTruth(iConfig)       ) ;
  childModules_.push_back(mod_trigger                          ) ;
  mod_trigger->config(this) ;
}

IIHEAnalysis::~IIHEAnalysis(){}

bool IIHEAnalysis::branchExists(std::string name){
  for(unsigned int i=0 ; i<allVars_.size() ; i++){
    if(allVars_.at(i)->name()==name) return true ;
  }
  return false ;
}
bool IIHEAnalysis::addBranch(std::string name){ return addBranch(name, currentVarType_) ; }

void IIHEAnalysis::setBranchType(int type){ currentVarType_ = type ; }
int  IIHEAnalysis::getBranchType(){ return currentVarType_ ; }

bool IIHEAnalysis::addBranch(std::string name, int type){
  // First check to see if this branch name has already been used
  bool success = !(branchExists(name)) ;
  if(success==false){
    return false ;
  }
  listOfBranches_.push_back(std::pair<std::string,int>(name,type)) ;
  switch(type){
    case kBool:{
      BranchWrapperB*  bw = new BranchWrapperB(name) ;
      vars_B_  .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kDouble:{
      BranchWrapperD*  bw = new BranchWrapperD(name) ;
      vars_D_  .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
    }
    case kFloat:{
      BranchWrapperF*  bw = new BranchWrapperF(name) ;
      vars_F_  .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kInt:{
      BranchWrapperI*  bw = new BranchWrapperI(name) ;
      vars_I_  .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kUInt:{
      BranchWrapperU*  bw = new BranchWrapperU(name) ;
      vars_U_  .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorBool:{
      BranchWrapperBV* bw = new BranchWrapperBV(name) ;
      vars_BV_ .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorDouble:{
      BranchWrapperDV* bw = new BranchWrapperDV(name) ;
      vars_DV_ .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorFloat:{
      BranchWrapperFV* bw = new BranchWrapperFV(name) ;
      vars_FV_ .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorInt:{
      BranchWrapperIV* bw = new BranchWrapperIV(name) ;
      vars_IV_ .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorUInt:{
      BranchWrapperUV* bw = new BranchWrapperUV(name) ;
      vars_UV_ .push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorBool:{
      BranchWrapperBVV* bw = new BranchWrapperBVV(name) ;
      vars_BVV_.push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorDouble:{
      BranchWrapperDVV* bw = new BranchWrapperDVV(name) ;
      vars_DVV_.push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorFloat:{
      BranchWrapperFVV* bw = new BranchWrapperFVV(name) ;
      vars_FVV_.push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorInt:{
      BranchWrapperIVV* bw = new BranchWrapperIVV(name) ;
      vars_IVV_.push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    case kVectorVectorUInt:{
      BranchWrapperUVV* bw = new BranchWrapperUVV(name) ;
      vars_UVV_.push_back(bw) ;
      allVars_.push_back((BranchWrapperBase*)bw) ;
      break ;
    }
    default : return false ; // Bail out if we don't know the type of branch
  }
  return true ;
}

// ------------ method called once each job just before starting event loop  -------------
void IIHEAnalysis::beginJob(){
  edm::Service<TFileService> fs;
  mainFile_ = new TFile("outfile.root", "RECREATE") ;
  dataTree_ = new TTree("IIHEAnalysis", "IIHEAnalysis") ;
  metaTree_ = new TTree("meta", "Information about globalTag etc") ;
  metaTree_->Branch("git_hash" , &git_hash_ ) ;
  metaTree_->Branch("globalTag", &globalTag_) ;
  metaTree_->Fill() ;
  for(unsigned int i=0 ; i<childModules_.size() ; i++){
    childModules_.at(i)->config(this) ;
    childModules_.at(i)->pubBeginJob() ;
  }
  configureBranches() ;
}

int IIHEAnalysis::saveToFile(TObject* obj){
  mainFile_->cd() ;
  return obj->Write() ;
}

void IIHEAnalysis::configureBranches(){
  for(unsigned int i=0 ; i<allVars_.size() ; i++){
    allVars_.at(i)->config(dataTree_) ;
  }
  return ;
}

void IIHEAnalysis::addToMCTruthWhitelist(std::vector<int> pdgIds){
  for(unsigned int i=0 ; i<pdgIds.size() ; i++){
    bool add = true ;
    for(unsigned int j=0 ; j<MCTruthWhitelist_.size() ; j++){
      if(abs(MCTruthWhitelist_.at(j))==abs(pdgIds.at(i))){
        add = false ;
        break ;
      }
    }
    if(add){
      MCTruthWhitelist_.push_back(pdgIds.at(i)) ;
    }
  }
}

bool IIHEAnalysis::addTriggerL1Electron(std::string name){
  for(unsigned int i=0 ; i<triggerL1FilterNamesElectron_.size() ; ++i){
    if(name==triggerL1FilterNamesElectron_.at(i)) return false ;
  }
  triggerL1FilterNamesElectron_.push_back(name) ;
  return true ;
}
bool IIHEAnalysis::addTriggerHLTElectron(std::string name, float DeltaRCut){
  for(unsigned int i=0 ; i<triggerHLTFilterNamesElectron_.size() ; ++i){
    if(name==triggerHLTFilterNamesElectron_.at(i).first) return false ;
  }
  std::pair<std::string,float> values(name,DeltaRCut) ;
  triggerHLTFilterNamesElectron_.push_back(values) ;
  return true ;
}

bool IIHEAnalysis::addTriggerL1Muon(std::string name){
  for(unsigned int i=0 ; i<triggerL1FilterNamesMuon_.size() ; ++i){
    if(name==triggerL1FilterNamesMuon_.at(i)) return false ;
  }
  triggerL1FilterNamesMuon_.push_back(name) ;
  return true ;
}
bool IIHEAnalysis::addTriggerHLTMuon(std::string name, float DeltaRCut){
  for(unsigned int i=0 ; i<triggerHLTFilterNamesMuon_.size() ; ++i){
    if(name==triggerHLTFilterNamesMuon_.at(i).first) return false ;
  }
  std::pair<std::string,float> values(name,DeltaRCut) ;
  triggerHLTFilterNamesMuon_.push_back(values) ;
  return true ;
}

// ------------ method called to for each event  -----------------------------------------

void IIHEAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  beginEvent() ;
  // Get the default collections
  // These should be harmonised across submodules, where possible
  iEvent.getByLabel(  photonCollectionLabel_,   photonCollection_) ;
  iEvent.getByLabel(electronCollectionLabel_, electronCollection_) ;
  iEvent.getByLabel(    muonCollectionLabel_,     muonCollection_) ;
  iEvent.getByLabel(     primaryVertexLabel_,       pvCollection_) ;
  
  // We take only the first primary vertex
  firstPrimaryVertex_ = new math::XYZPoint(0.0,0.0,0.0) ;
  const reco::VertexCollection* primaryVertices = getPrimaryVertices() ;
  if(primaryVertices->size()>0){
    reco::VertexCollection::const_iterator firstpv = primaryVertices->begin();
    firstPrimaryVertex_->SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  for(unsigned i=0 ; i<childModules_.size() ; i++){
    childModules_.at(i)->pubAnalyze(iEvent, iSetup) ;
  }
  endEvent() ;
  dataTree_->Fill() ;
}

void IIHEAnalysis::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}

void IIHEAnalysis::beginEvent(){
  for(unsigned int i=0 ; i<childModules_.size() ; i++){ childModules_.at(i)->pubBeginEvent() ; }
  for(unsigned int i=0 ; i<allVars_.size()      ; i++){ allVars_.at(i)->beginEvent()         ; }
}

void IIHEAnalysis::endEvent(){
  for(unsigned int i=0 ; i<childModules_.size() ; i++){ childModules_.at(i)->pubEndEvent() ; }
  for(unsigned int i=0 ; i<allVars_.size()      ; i++){      allVars_.at(i)->endEvent()    ; }
}


// ------------ method called once each job just after ending the event loop  ------------
void 
IIHEAnalysis::endJob(){
  std::vector<std::string> untouchedBranchNames ;
  for(unsigned int i=0 ; i<allVars_.size() ; i++){
    if(allVars_.at(i)->is_touched()==false) untouchedBranchNames.push_back(allVars_.at(i)->name()) ;
  }
  if(debug_==true){
    if(untouchedBranchNames.size()>0){
      std::cout << "The following branches were never touched:" << std::endl ;
      for(unsigned int i=0 ; i<untouchedBranchNames.size() ; i++){
        std::cout << "  " << untouchedBranchNames.at(i) << std::endl ;
      }
    }
  }
  
  if(mainFile_){
    mainFile_->Write() ;
    delete mainFile_ ;
  }
}

bool IIHEAnalysis::store(std::string name, bool value){
  for(unsigned int i=0 ; i<vars_B_.size() ; i++){
    if(vars_B_.at(i)->name()==name){
      vars_B_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_BV_.size() ; i++){
    if(vars_BV_.at(i)->name()==name){
      vars_BV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (bool) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, double value){
  // Try to fill doubles, then floats
  for(unsigned int i=0 ; i<vars_D_.size() ; i++){
    if(vars_D_.at(i)->name()==name){
      vars_D_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){
    if(vars_DV_.at(i)->name()==name){
      vars_DV_ .at(i)->push(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_F_.size() ; i++){
    if(vars_F_.at(i)->name()==name){
      vars_F_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){
    if(vars_FV_.at(i)->name()==name){
      vars_FV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (double) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, float value){
  // Try to fill floats, then doubles
  for(unsigned int i=0 ; i<vars_F_.size() ; i++){
    if(vars_F_.at(i)->name()==name){
      vars_F_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){
    if(vars_FV_.at(i)->name()==name){
      vars_FV_ .at(i)->push(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_D_.size() ; i++){
    if(vars_D_.at(i)->name()==name){
      vars_D_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){
    if(vars_DV_.at(i)->name()==name){
      vars_DV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (float) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, int value){
  for(unsigned int i=0 ; i<vars_I_.size() ; i++){
    if(vars_I_.at(i)->name()==name){
      vars_I_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){
    if(vars_IV_.at(i)->name()==name){
      vars_IV_ .at(i)->push(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_U_.size() ; i++){
    if(vars_U_.at(i)->name()==name){
      vars_U_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_UV_.size() ; i++){
    if(vars_UV_.at(i)->name()==name){
      vars_UV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (int) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, unsigned int value){
  for(unsigned int i=0 ; i<vars_U_.size() ; i++){
    if(vars_U_.at(i)->name()==name){
      vars_U_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_UV_.size() ; i++){
    if(vars_UV_.at(i)->name()==name){
      vars_UV_ .at(i)->push(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_I_.size() ; i++){
    if(vars_I_.at(i)->name()==name){
      vars_I_ .at(i)->set(value) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){
    if(vars_IV_.at(i)->name()==name){
      vars_IV_ .at(i)->push(value) ;
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (uint) branch named " << name << std::endl ;
  return false ;
}

bool IIHEAnalysis::store(std::string name, std::vector<bool> values){
  for(unsigned int i=0 ; i<vars_BVV_.size() ; i++){
    if(vars_BVV_.at(i)->name()==name){
      vars_BVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_BV_.size() ; i++){
    if(vars_BV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_BV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (vector bool) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, std::vector<float> values){
  for(unsigned int i=0 ; i<vars_FVV_.size() ; i++){
    if(vars_FVV_.at(i)->name()==name){
      vars_FVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_FV_.size() ; i++){
    if(vars_FV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_FV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (vector float) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, std::vector<double> values){
  for(unsigned int i=0 ; i<vars_DVV_.size() ; i++){
    if(vars_DVV_.at(i)->name()==name){
      vars_DVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_DV_.size() ; i++){
    if(vars_DV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_DV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (vector double) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, std::vector<int> values){
  for(unsigned int i=0 ; i<vars_IVV_.size() ; i++){
    if(vars_IVV_.at(i)->name()==name){
      vars_IVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_IV_.size() ; i++){
    if(vars_IV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_IV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (vector int) branch named " << name << std::endl ;
  return false ;
}
bool IIHEAnalysis::store(std::string name, std::vector<unsigned int> values){
  for(unsigned int i=0 ; i<vars_UVV_.size() ; i++){
    if(vars_UVV_.at(i)->name()==name){
      vars_UVV_ .at(i)->push(values) ;
      return true ;
    }
  }
  for(unsigned int i=0 ; i<vars_UV_.size() ; i++){
    if(vars_UV_.at(i)->name()==name){
      for(unsigned j=0 ; j<values.size() ; j++){
        vars_UV_ .at(i)->push(values.at(j)) ;
      }
      return true ;
    }
  }
  if(debug_) std::cout << "Could not find a (vector uint) branch named " << name << std::endl ;
  return false ;
}

DEFINE_FWK_MODULE(IIHEAnalysis);
