#ifndef UserCode_IIHETree_TriggerObject_h
#define UserCode_IIHETree_TriggerObject_h

#include "UserCode/IIHETree/interface/Types.h"
#include "UserCode/IIHETree/interface/IIHEAnalysis.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Framework/interface/Event.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std ;
using namespace reco;
using namespace edm ;

class IIHEAnalysis ;

class TriggerMatchParameters{
private:
  int triggerLevel_ ;
  int particleType_ ;
  std::string triggerName_ ; // Name of the trigger
  std::string prefix_      ; // Used to identify the analysis that wants this matching (eg HEEP_)
  std::string filterName_  ; // Name of the EDAnalyzer etc that made the trigger decision
  std::string branchName_  ; // Name of the branch, which must be unique
  int filterIndex_         ; // Index of the filter used for DeltaR matching
public:
  TriggerMatchParameters(int, int, std::string, std::string) ;
  TriggerMatchParameters(TriggerMatchParameters*, std::string) ;
  ~TriggerMatchParameters(){} ;
  TriggerMatchParameters* Clone() ;
  std::string      prefix(){ return prefix_      ; }
  std::string triggerName(){ return triggerName_ ; }
  std::string  filterName(){ return filterName_  ; }
  std::string  branchName(){ return "trigMatch_" + prefix_ + triggerName_ + "_" + filterName_ + "_DeltaR" ; }
  const int particleType(){ return particleType_ ; }
  const int triggerLevel(){ return triggerLevel_ ; }
  int setFilterIndex(edm::Handle<trigger::TriggerEvent>, edm::InputTag) ;
  void setFilterName(std::string filterName){ filterName_ = filterName ; }
  float matchObject(edm::Handle<trigger::TriggerEvent>, float, float) ;
};

class L1Trigger{
private:
  std::string name_ ;
  std::string branchName_ ;
  int filterIndex_;
  bool accept_  ;
  bool touched_ ;
  int prescale_ ;
  int index_ ;
  
  double barrelEnd_       ;
  double regionEtaSizeEB_ ;
  double regionEtaSizeEE_ ;
  double regionPhiSize_   ;
  
  bool matchObject(edm::Handle<trigger::TriggerEvent>, float, float) ;
public:
  L1Trigger(std::string, std::string) ;
  ~L1Trigger() ;
  void reset() ;
  
  std::string       name(){ return       name_ ; }
  std::string branchName(){ return branchName_ ; }
  int setFilterIndex(edm::Handle<trigger::TriggerEvent>, edm::InputTag) ;
  bool matchElectron(edm::Handle<trigger::TriggerEvent>, reco::GsfElectron*) ;
  bool matchMuon    (edm::Handle<trigger::TriggerEvent>, reco::Muon*       ) ;
};


class HLTrigger{
private:
  std::string name_ ;
  int  accept_  ;
  bool touched_ ;
  bool error_ ;
  int  prescale_ ;
  int  index_ ;
  int  searchStatus_ ;
  
  int nSC_    ;
  int nPh_    ;
  int nEl_    ;
  int nMu_    ;
  int nTau_   ;
  int nJet_   ;
  int hasMET_ ;
  int nSCEl_  ;
  int nSCPh_  ;
  int nTypes_ ;
  
  std::string acceptBranchName_ ;
  std::string prescaleBranchName_ ;
  
  enum searchStatuses{ notSearchedFor , searchedForAndFound , searchedForAndNotFound } ;
  std::vector<TriggerMatchParameters*> matchingParameters_ ;
  
  int nSubstringInString(const std::string&, const std::string&) ;
  int nMuonsInTriggerName() ;
  int nSuperclustersInTriggerName() ;
  int nPhotonsInTriggerName() ;
  int nElectronsInTriggerName() ;
  int nTausInTriggerName() ;
  int nJetsInTriggerName() ;
  int METInTriggerName() ;
  
public:
  HLTrigger(std::string) ;
  ~HLTrigger() ;
  void reset() ;
  int createBranches(IIHEAnalysis*, int) ;
  bool findIndex(std::vector<std::string>) ;
  bool  beginRun(std::vector<std::string>) ;
  
  bool status(const edm::Event&, edm::EventSetup const&, HLTConfigProvider const&, Handle<TriggerResults> const&, edm::Handle<trigger::TriggerEvent> const&, edm::InputTag const&) ;
  void fill(IIHEAnalysis*, edm::Handle<trigger::TriggerEvent>, edm::InputTag) ;
  bool addMatching(TriggerMatchParameters*) ;
  
  std::string name(){ return name_ ; }
  void setIndex(int index){ index_ = index ; }
  int index(){ return index_ ; }
  void printFilterNames() ;
  
  bool isSingleElectron(){ return nEl_==1 ; }
  bool isDoubleElectron(){ return nEl_==2 ; }
  bool isTripleElectron(){ return nEl_==3 ; }
  bool isSingleMuon(){ return nMu_==1 ; }
  bool isDoubleMuon(){ return nMu_==2 ; }
  bool isTripleMuon(){ return nMu_==3 ; }
  bool isSingleElectronSingleMuon(){ return (nEl_==1 && nMu_==1) ; }
  bool isSingleElectronDoubleMuon(){ return (nEl_==1 && nMu_==2) ; }
  bool isDoubleElectronSingleMuon(){ return (nEl_==2 && nMu_==1) ; }
  bool isOnlySingleElectron(){ return (nTypes_==1*pow(10,(int)kElectron)) ; }
  bool isOnlyDoubleElectron(){ return (nTypes_==1*pow(10,(int)kElectron)) ; }
  bool isOnlyTripleElectron(){ return (nTypes_==1*pow(10,(int)kElectron)) ; }
  bool isOnlySingleMuon(){ return (nTypes_==1*pow(10,(int)kMuon)) ; }
  bool isOnlyDoubleMuon(){ return (nTypes_==1*pow(10,(int)kMuon)) ; }
  bool isOnlyTripleMuon(){ return (nTypes_==1*pow(10,(int)kMuon)) ; }
  bool isOnlySingleElectronSingleMuon(){ return (nTypes_ = 1*pow(10,(int)kElectron) + 1*pow(10,(int)kMuon)) ; }
  bool isOnlySingleElectronDoubleMuon(){ return (nTypes_ = 1*pow(10,(int)kElectron) + 2*pow(10,(int)kMuon)) ; }
  bool isOnlyDoubleElectronSingleMuon(){ return (nTypes_ = 2*pow(10,(int)kElectron) + 1*pow(10,(int)kMuon)) ; }
  
};

#endif
