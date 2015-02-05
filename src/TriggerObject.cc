#include "UserCode/IIHETree/interface/TriggerObject.h"

TriggerMatchParameters::TriggerMatchParameters(int triggerLevel, int particleType, std::string triggerName, std::string filterName){
  triggerLevel_ = triggerLevel ;
  particleType_ = particleType ;
  triggerName_  = triggerName ;
  filterName_   = filterName ;
  prefix_ = "" ;
  switch(particleType_){
    case kSuperCluster: prefix_ = "SC_"  ; break ;
    case kPhoton      : prefix_ = "ph_"  ; break ;
    case kElectron    : prefix_ = "el_"  ; break ;
    case kMuon        : prefix_ = "mu_"  ; break ;
    case kTau         : prefix_ = "tau_" ; break ;
    case kJet         : prefix_ = "jet_" ; break ;
  }
  filterIndex_ = -1 ;
}
TriggerMatchParameters::TriggerMatchParameters(TriggerMatchParameters* TMP, std::string filterName){
  triggerLevel_ = TMP->triggerLevel() ;
  particleType_ = TMP->particleType() ;
  triggerName_  = TMP->triggerName()  ;
  filterName_   = filterName ;
  prefix_ = "" ;
  switch(particleType_){
    case kSuperCluster: prefix_ = "SC_"  ; break ;
    case kPhoton      : prefix_ = "ph_"  ; break ;
    case kElectron    : prefix_ = "el_"  ; break ;
    case kMuon        : prefix_ = "mu_"  ; break ;
    case kTau         : prefix_ = "tau_" ; break ;
    case kJet         : prefix_ = "jet_" ; break ;
  }
  branchName_   = "trigMatch_" + prefix_ + triggerName_ + "_" + filterName_ + "_DeltaR" ;
  std::cout << prefix_ << " " << particleType_ << " " << branchName_ << std::endl ;
  filterIndex_ = -1 ;
}

TriggerMatchParameters* TriggerMatchParameters::Clone(){
  TriggerMatchParameters* TMP = new TriggerMatchParameters(triggerLevel_, particleType_, triggerName_, filterName_) ;
  return TMP ;
}


int TriggerMatchParameters::setFilterIndex(edm::Handle<trigger::TriggerEvent> trigEvent, edm::InputTag trigEventTag){
  filterIndex_ = trigEvent->filterIndex(edm::InputTag(filterName_,"",trigEventTag.process())) ;
  return filterIndex_ ;
}
float TriggerMatchParameters::matchObject(edm::Handle<trigger::TriggerEvent> trigEvent, float eta, float phi){
  if(filterIndex_<0) return 20 ;
  
  if(filterIndex_<trigEvent->sizeFilters()){ 
    const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex_) ; 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects()) ;
    
    // Now loop over the trigger objects passing filter
    float smallestDR = 999 ;
    for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt!=trigKeys.end(); ++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt] ;
      float DR = deltaR(eta,phi,obj.eta(),obj.phi()) ;
      if(DR<smallestDR) smallestDR = DR ;
    }
    return smallestDR ;
  }
  return 10 ;
}

L1Trigger::L1Trigger(std::string name, std::string prefix){
  filterIndex_ = -1 ;
  barrelEnd_ = 1.4791 ;
  regionEtaSizeEB_ = 0.522 ;
  regionEtaSizeEE_ = 1.0 ;
  regionPhiSize_ = 1.044 ;
  name_ = name ;
  branchName_ = prefix + name_ ;
  reset() ;
}
L1Trigger::~L1Trigger(){}
void L1Trigger::reset(){
  accept_ = false ;
  touched_ = false ;
  prescale_ = -999 ;
}

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

HLTrigger::HLTrigger(std::string name){
  name_ = name ;
  index_ = -1 ;
  searchStatus_ = notSearchedFor ;
  reset() ;
  acceptBranchName_   = "trig_" + name + "_accept"   ;
  prescaleBranchName_ = "trig_" + name + "_prescale" ;
  nSC_    = nSuperclustersInTriggerName() ;
  nPh_    = nPhotonsInTriggerName() ;
  nEl_    = nElectronsInTriggerName() ;
  nMu_    = nMuonsInTriggerName() ;
  nTau_   = nTausInTriggerName() ;
  nJet_   = nJetsInTriggerName() ;
  hasMET_ = METInTriggerName() ;
  nSCEl_  = nSC_+nEl_ ;
  
  // Slightly easier way to handle multiple objects
  nTypes_ =    nSC_*pow(10,(int)kSuperCluster)
          +    nPh_*pow(10,(int)kPhoton)
          +    nEl_*pow(10,(int)kElectron)
          +    nMu_*pow(10,(int)kMuon)
          +   nTau_*pow(10,(int)kTau)
          +   nJet_*pow(10,(int)kJet)
          + hasMET_*pow(10,(int)kMET) ;
}
HLTrigger::~HLTrigger(){}
void HLTrigger::reset(){
  accept_ = -1 ;
  touched_ = false ;
  prescale_ = -999 ;
}

int HLTrigger::nSuperclustersInTriggerName(){
  int scCount = nSubstringInString(name_, "_SC") ;
  return scCount ;
}
int HLTrigger::nPhotonsInTriggerName(){
  int scCount = nSubstringInString(name_, "Ph") ;
  return scCount ;
}
int HLTrigger::nElectronsInTriggerName(){
  int singleElectronCount = nSubstringInString(name_, "Ele"      ) ;
  int doubleElectronCount = nSubstringInString(name_, "DoubleEle") ;
  int tripleElectronCount = nSubstringInString(name_, "TripleEle") ;
  int totalElectronCount = 3*tripleElectronCount + 2*doubleElectronCount + singleElectronCount ;
  return totalElectronCount ;
}

int HLTrigger::nMuonsInTriggerName(){
  int singleMuonCount = nSubstringInString(name_, "Mu"      ) ;
  int doubleMuonCount = nSubstringInString(name_, "DoubleMu") ;
  int tripleMuonCount = nSubstringInString(name_, "TripleMu") ;
  int totalMuonCount = 3*tripleMuonCount + 2*doubleMuonCount + singleMuonCount ;
  return totalMuonCount ;
}
int HLTrigger::nTausInTriggerName(){
  int tauCount = nSubstringInString(name_, "Tau") ;
  return tauCount ;
}
int HLTrigger::nJetsInTriggerName(){
  int jetCount = nSubstringInString(name_, "Jet") ;
  return jetCount ;
}
int HLTrigger::METInTriggerName(){
  int metCount = nSubstringInString(name_, "MET") ;
  return metCount ;
}
int HLTrigger::nSubstringInString(const std::string& str, const std::string& sub){
  // Taken from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/HLTriggerOffline/Egamma/src/EmDQM.cc#L1064
  // Thanks, Thomas!
  if(sub.length()==0) return 0 ;
  int count = 0 ;
  for (size_t offset=str.find(sub) ; offset!=std::string::npos ; offset=str.find(sub, offset + sub.length())){ ++count ; }
  return count;
}


bool HLTrigger::status(const edm::Event& iEvent, edm::EventSetup const& iSetup, HLTConfigProvider const& hltConfig, Handle<TriggerResults> const& triggerResults, edm::Handle<trigger::TriggerEvent> const& trigEvent, edm::InputTag const& trigEventTag){
  if(searchStatus_==searchedForAndFound && index_>=0){
    touched_  = true ;
    accept_   = triggerResults->accept(index_) ;
    prescale_ = hltConfig.prescaleValue(iEvent, iSetup, name_) ;
    return true ;
  }
  else{
    return false ;
  }
}
void HLTrigger::fill(IIHEAnalysis* analysis, edm::Handle<trigger::TriggerEvent> trigEvent, edm::InputTag trigEventTag){
  analysis->store(acceptBranchName_  , accept_  ) ;
  analysis->store(prescaleBranchName_, prescale_) ;
  
  for(unsigned int i=0 ; i<matchingParameters_.size() ; ++i){
    TriggerMatchParameters* MP = matchingParameters_.at(i) ;
    MP->setFilterIndex(trigEvent, trigEventTag) ;
    std::string matchBranchName = MP->branchName() ;
    
    switch(MP->particleType()){
      case kSuperCluster:{
        std::vector<const reco::SuperCluster*> superclusters = analysis->getSuperClusters() ;
        for(unsigned int i_sc=0 ; i_sc<superclusters.size() ; i_sc++){
          reco::SuperCluster* sc = (reco::SuperCluster*)superclusters.at(i_sc) ;
          float DeltaR = MP->matchObject(trigEvent, sc->eta(), sc->phi()) ;
          analysis->store(matchBranchName, DeltaR) ;
        }
      }
      case kPhoton:{
        reco::PhotonCollection photons = analysis->getPhotonCollection() ;
        for(reco::PhotonCollection::const_iterator phiter = photons.begin() ; phiter!=photons.end() ; ++phiter){
          float DeltaR = MP->matchObject(trigEvent, phiter->eta(), phiter->phi()) ;
          analysis->store(matchBranchName, DeltaR) ;
        }
      }
      case kElectron:{
        reco::GsfElectronCollection electrons = analysis->getElectronCollection() ;
        for(reco::GsfElectronCollection::const_iterator gsfiter=electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
          float DeltaR = MP->matchObject(trigEvent, gsfiter->eta(), gsfiter->phi()) ;
          analysis->store(matchBranchName, DeltaR) ;
        }
      }
      case kMuon:{
        reco::MuonCollection muons = analysis->getMuonCollection() ;
        for(reco::MuonCollection::const_iterator muIt = muons.begin(); muIt != muons.end(); ++muIt){
          float DeltaR = MP->matchObject(trigEvent, muIt->eta(), muIt->phi()) ;
          analysis->store(matchBranchName, DeltaR) ;
        }
      }
    }
  }
}

int HLTrigger::createBranches(IIHEAnalysis* analysis, int nEvents){
  int result = 0 ;
  bool result_acceptBranch   = analysis->addBranch(  acceptBranchName_, kInt) ;
  bool result_prescaleBranch = analysis->addBranch(prescaleBranchName_, kInt) ;
  
  if(result_acceptBranch){
    result++ ;
    for(int i=0 ; i<nEvents ; ++i){
      analysis->store(acceptBranchName_, -1) ;
    }
  }
  if(result_prescaleBranch){
    result++ ;
    for(int i=0 ; i<nEvents ; ++i){
      analysis->store(prescaleBranchName_, -1) ;
    }
  }
  
  for(unsigned int i=0 ; i<matchingParameters_.size() ; ++i){
    TriggerMatchParameters* MP  = matchingParameters_.at(i) ;
    std::string matchBranchName = MP->branchName() ;
    std::cout << matchBranchName << std::endl ;
    
    bool result_matchBranch = analysis->addBranch(matchBranchName, kVectorFloat) ;
    if(result_matchBranch){
      result++ ;
      for(int i=0 ; i<nEvents ; ++i){
        analysis->store(matchBranchName, -1) ;
      }
    }
  }
  
  return result ;
}

bool HLTrigger::beginRun(std::vector<std::string> names){
  bool success = findIndex(names) ;
  return success ;
}
bool HLTrigger::findIndex(std::vector<std::string> names){
  searchStatus_ = notSearchedFor ;
  for(unsigned int i=0 ; i<names.size() ; ++i){
    if(names.at(i)==name_){
      index_ = i ;
      searchStatus_ = searchedForAndFound ;
      return true ;
    }
  }
  index_ = -1 ;
  searchStatus_ = searchedForAndNotFound ;
  return false ;
}
bool HLTrigger::addMatching(TriggerMatchParameters* TMP){
  if(TMP->triggerName()!=name_) return false ;
  for(unsigned int i=0 ; i<matchingParameters_.size() ; ++i){
    TriggerMatchParameters* MP = matchingParameters_.at(i) ;
    if(TMP->branchName()==MP->branchName()){
      return false ;
    }
  }
  matchingParameters_.push_back(TMP) ;
  return true ;
}
void HLTrigger::printFilterNames(){
  //for(unsigned int i=0 ; i<moduleNames.size() ; ++i){
  //  std::cout << i << " " << moduleNames.at(i) << std::endl ;
  //}
}


