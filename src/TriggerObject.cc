#include "UserCode/IIHETree/interface/TriggerObject.h"

TriggerFilter::TriggerFilter(std::string name, std::string triggerName){
    name_ = name ;
    triggerName_ = triggerName ;
    etaBranchName_ = "trig_" + triggerName_ + "_" + name_ + "_eta" ;
    phiBranchName_ = "trig_" + triggerName_ + "_" + name_ + "_phi" ;
}
int TriggerFilter::createBranches(IIHEAnalysis* analysis){
  int result = 0 ;
  result += analysis->addBranch(etaBranchName_, kVectorFloat) ;
  result += analysis->addBranch(phiBranchName_, kVectorFloat) ;
  return result ;
}
int TriggerFilter::setIndex(edm::Handle<trigger::TriggerEvent> trigEvent, edm::InputTag trigEventTag){
  index_ = trigEvent->filterIndex(edm::InputTag(name_,"",trigEventTag.process())) ;
  return index_ ;
}
int TriggerFilter::setValues(edm::Handle<trigger::TriggerEvent> trigEvent, IIHEAnalysis* analysis){
  etaValues_.clear() ;
  phiValues_.clear() ;
  
  if(index_<0) return 2 ;
  if(index_<trigEvent->sizeFilters()){
    const trigger::Keys& trigKeys = trigEvent->filterKeys(index_) ; 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects()) ;
    
    // Now loop over the trigger objects passing filter
    for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt!=trigKeys.end(); ++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt] ;
      analysis->store(etaBranchName_, obj.eta()) ;
      analysis->store(phiBranchName_, obj.phi()) ;
      etaValues_.push_back(obj.eta()) ;
      phiValues_.push_back(obj.phi()) ;
    }
    return 0 ;
  }
  return 1 ;
}
bool TriggerFilter::store(IIHEAnalysis* analysis){
  bool etaSuccess = analysis->store(etaBranchName_, etaValues_) ;
  bool phiSuccess = analysis->store(phiBranchName_, phiValues_) ;
  return (etaSuccess && phiSuccess) ;
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

HLTrigger::HLTrigger(std::string name, HLTConfigProvider hltConfig){
  name_ = name ;
  index_ = -1 ;
  searchStatus_ = notSearchedFor ;
  reset() ;
  acceptBranchName_   = "trig_" + name + "_accept"   ;
  prescaleBranchName_ = "trig_" + name + "_prescale" ;
  nSC_     = nSuperclustersInTriggerName() ;
  nPh_     = nPhotonsInTriggerName() ;
  nEl_     = nElectronsInTriggerName() ;
  nMu_     = nMuonsInTriggerName() ;
  nTau_    = nTausInTriggerName() ;
  nJet_    = nJetsInTriggerName() ;
  nBJet_   = nBJetsInTriggerName() ;
  hasMET_  = METInTriggerName() ;
  hasHT_   = HTInTriggerName() ;
  hasALCa_ = ALCaInTriggerName() ;
  nSCEl_   = nSC_+nEl_ ;
  
  // Slightly easier way to handle multiple objects
  nTypes_ =     nSC_*pow(10,(int)kSuperCluster)
          +     nPh_*pow(10,(int)kPhoton)
          +     nEl_*pow(10,(int)kElectron)
          +     nMu_*pow(10,(int)kMuon)
          +    nTau_*pow(10,(int)kTau)
          +    nJet_*pow(10,(int)kJet)
          +   nBJet_*pow(10,(int)kBJet)
          +  hasMET_*pow(10,(int)kMET) 
          +   hasHT_*pow(10,(int)kHT) 
          + hasALCa_*pow(10,(int)kALCa) ;
  
  findIndex(hltConfig) ;
  std::vector<std::string> moduleNames = hltConfig.moduleLabels(index_) ;
  std::vector<std::string> moduleNamesWithTags ;
  for(unsigned int j=0 ; j<moduleNames.size() ; ++j){
    if(hltConfig.saveTags(moduleNames.at(j))){
      moduleNamesWithTags.push_back(moduleNames.at(j)) ;
    }
  }
  for(unsigned int i=0 ; i<moduleNamesWithTags.size() ; ++i){
    filters_.push_back(new TriggerFilter(moduleNamesWithTags.at(i), name_)) ;
  }
  
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
  int singlePhotonCount = nSubstringInString(name_, "Photon"      ) ;
  int doublePhotonCount = nSubstringInString(name_, "DoublePhoton") ;
  int triplePhotonCount = nSubstringInString(name_, "TriplePhoton") ;
  int totalPhotonCount = 2*triplePhotonCount + 1*doublePhotonCount + singlePhotonCount ;
  return totalPhotonCount ;
}
int HLTrigger::nElectronsInTriggerName(){
  int singleElectronCount = nSubstringInString(name_, "Ele"      ) ;
  int doubleElectronCount = nSubstringInString(name_, "DoubleEle") + nSubstringInString(name_, "DiEle") ;
  int tripleElectronCount = nSubstringInString(name_, "TripleEle") ;
  int totalElectronCount = 2*tripleElectronCount + 1*doubleElectronCount + singleElectronCount ;
  return totalElectronCount ;
}

int HLTrigger::nMuonsInTriggerName(){
  int singleMuonCount = nSubstringInString(name_, "Mu"      ) ;
  int doubleMuonCount = nSubstringInString(name_, "DoubleMu") + nSubstringInString(name_, "DiMu") + nSubstringInString(name_, "Dimuon") ;
  int tripleMuonCount = nSubstringInString(name_, "TripleMu") ;
  int totalMuonCount = 2*tripleMuonCount + 1*doubleMuonCount + singleMuonCount ;
  return totalMuonCount ;
}
int HLTrigger::nTausInTriggerName(){
  int tauCount = nSubstringInString(name_, "Tau") ;
  return tauCount ;
}
int HLTrigger::nJetsInTriggerName(){
  int ignoreJetCount = nSubstringInString(name_, "NoJetId"  ) ;
  int singleJetCount = nSubstringInString(name_, "Jet"  ) ;
  int doubleJetCount = nSubstringInString(name_, "DiJet") + nSubstringInString(name_, "DiPFJet") + nSubstringInString(name_, "DoubleJet") + nSubstringInString(name_, "DiCentralJet") + nSubstringInString(name_, "DiCentralPFJet") ;
  int tripleJetCount = nSubstringInString(name_, "TriCentralPFJet" ) ;
  int quadJetCount  = nSubstringInString(name_, "QuadJet" ) + nSubstringInString(name_, "QuadPFJet" ) ;
  int sixJetCount   = nSubstringInString(name_, "SixJet"  ) + nSubstringInString(name_, "SixPFJet"  ) ;
  int eightJetCount = nSubstringInString(name_, "EightJet") + nSubstringInString(name_, "EightPFJet") ;
  int totalJetCount = 7*eightJetCount + 5*sixJetCount + 3*quadJetCount + 2*tripleJetCount + 1*doubleJetCount + singleJetCount - 1*ignoreJetCount ;
  return totalJetCount ;
}
int HLTrigger::nBJetsInTriggerName(){
  int singleBJetCount = nSubstringInString(name_, "BJet"  ) ;
  int doubleBJetCount = nSubstringInString(name_, "DiBJet"  ) + nSubstringInString(name_,  "DiCentral") ;
  int tripleBJetCount = nSubstringInString(name_, "TriiBJet") + nSubstringInString(name_, "TriCentral") ;
  int totalBJetCount = 2*tripleBJetCount + 1*doubleBJetCount + singleBJetCount ;
  return totalBJetCount ;
}
int HLTrigger::METInTriggerName(){
  int METCount = nSubstringInString(name_, "MET") ;
  return METCount ;
}
int HLTrigger::HTInTriggerName(){
  int HTCount = nSubstringInString(name_, "HT") ;
  return HTCount ;
}
int HLTrigger::ALCaInTriggerName(){
  int ALCaCount = nSubstringInString(name_, "ALCa") ;
  return ALCaCount ;
}
int HLTrigger::nSubstringInString(const std::string& str, const std::string& sub){
  // Taken from https://github.com/cms-sw/cmssw/blob/CMSSW_7_4_X/HLTriggerOffline/Egamma/src/EmDQM.cc#L1064
  // Thanks, Thomas!
  if(sub.length()==0) return 0 ;
  int count = 0 ;
  for (size_t offset=str.find(sub) ; offset!=std::string::npos ; offset=str.find(sub, offset + sub.length())){ ++count ; }
  return count;
}

int HLTrigger::status(const edm::Event& iEvent, edm::EventSetup const& iSetup, HLTConfigProvider const& hltConfig, Handle<TriggerResults> const& triggerResults, edm::Handle<trigger::TriggerEvent> trigEvent, IIHEAnalysis* analysis){
  if(searchStatus_==searchedForAndFound && index_>=0){
    touched_  = true ;
    accept_   = triggerResults->accept(index_) ;
    prescale_ = hltConfig.prescaleValue(iEvent, iSetup, name_) ;
    for(unsigned i=0 ; i<filters_.size() ; ++i){
      filters_.at(i)->setValues(trigEvent, analysis) ;
    }
    return 0 ;
  }
  return 2 ;
}
void HLTrigger::store(IIHEAnalysis* analysis){
  analysis->store(  acceptBranchName_, accept_  ) ;
  analysis->store(prescaleBranchName_, prescale_) ;
}

int HLTrigger::createBranches(IIHEAnalysis* analysis){
  int result = 0 ;
  result += analysis->addBranch(  acceptBranchName_, kInt) ;
  result += analysis->addBranch(prescaleBranchName_, kInt) ;
  
  for(unsigned i=0 ; i<filters_.size() ; ++i){
    result += filters_.at(i)->createBranches(analysis) ;
  }
  
  return result ;
}

bool HLTrigger::beginRun(HLTConfigProvider const& hltConfig){
  bool success = findIndex(hltConfig) ;
  return success ;
}
int HLTrigger::findIndex(HLTConfigProvider const& hltConfig){
  searchStatus_ = notSearchedFor ;
  std::vector<std::string> names = hltConfig.triggerNames() ;
  for(unsigned int i=0 ; i<names.size() ; ++i){
    if(names.at(i)==name_){
      index_ = i ;
      searchStatus_ = searchedForAndFound ;
      return 0 ;
    }
  }
  index_ = -1 ;
  searchStatus_ = searchedForAndNotFound ;
  return 1 ;
}
bool HLTrigger::addFilter(std::string fileName){
  return true ;
}


