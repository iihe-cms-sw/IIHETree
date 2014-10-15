#ifndef UserCode_IIHETree_TriggerObject_h
#define UserCode_IIHETree_TriggerObject_h

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std ;
using namespace reco;
using namespace edm ;

class L1Trigger{
private:
  std::string name_ ;
  std::string branchName_ ;
  int filterIndex_;
  
  const double barrelEnd_       ;
  const double regionEtaSizeEB_ ;
  const double regionEtaSizeEE_ ;
  const double regionPhiSize_   ;
  
  bool matchObject(edm::Handle<trigger::TriggerEvent>, float, float) ;
public:
  L1Trigger(std::string, std::string) ;
  ~L1Trigger() ;
  
  std::string       name(){ return       name_ ; }
  std::string branchName(){ return branchName_ ; }
  int setFilterIndex(edm::Handle<trigger::TriggerEvent>, edm::InputTag) ;
  bool matchElectron(edm::Handle<trigger::TriggerEvent>, reco::GsfElectron*) ;
  bool matchMuon    (edm::Handle<trigger::TriggerEvent>, reco::Muon*       ) ;
};

class HLTrigger{
private:
  std::string name_ ;
  std::string branchName_ ;
  int filterIndex_;
  float DeltaRCut_ ;
  
  bool matchObject(edm::Handle<trigger::TriggerEvent>, float, float) ;
public:
  HLTrigger(std::string, std::string, float) ;
  ~HLTrigger() ;
  
  std::string       name(){ return       name_ ; }
  std::string branchName(){ return branchName_ ; }
  int setFilterIndex(edm::Handle<trigger::TriggerEvent>, edm::InputTag) ;
  bool matchElectron(edm::Handle<trigger::TriggerEvent>, reco::GsfElectron*) ;
  bool matchMuon    (edm::Handle<trigger::TriggerEvent>, reco::Muon*       ) ;
};

#endif
