#ifndef UserCode_IIHETree_HEEPCut_h
#define UserCode_IIHETree_HEEPCut_h

#include <algorithm>

#include "UserCode/IIHETree/interface/IIHEModuleHEEP.h"
// These classes are used to make it easier to keep track of HEEP cuts
// Each class corresponds to a different cut
// The cuts can be arranged into cut collections (eg HEEP41, HEEP50) to process a whole list of cuts at once

class IIHEModuleHEEP ;

class HEEPCutBase{
private:
  std::string name_ ;
  std::string branchName_n_ ;
  std::string branchName_nCumulative_ ;
  std::string branchName_value_ ;
  bool status_ ;
  int index_ ;
  int nPass_ ;
  int nPassCumulative_ ;
  IIHEModuleHEEP* parent_mod_ ;
public:
  HEEPCutBase(std::string, IIHEModuleHEEP*) ;
  ~HEEPCutBase() ;
  virtual bool applyCut(reco::GsfElectron*, bool) ;
  bool addBranches() ;
  void store() ;
  void reset() ;
  void setStatus(bool, bool) ;
  bool getStatus() ;
  std::string name(){ return name_ ; }
  
  float value_ ; // Value of the variable.  public for now, should be changed to private with accessor methods later
  
  void beginEvent() ;
  void endEvent() ;
  
  bool isBarrel(reco::GsfElectron*) ;
  bool isEndcap(reco::GsfElectron*) ;
  int detectorRegion(reco::GsfElectron*) ;
  
  bool isBarrel_41(reco::GsfElectron*) ;
  bool isEndcap_41(reco::GsfElectron*) ;
  int detectorRegion_41(reco::GsfElectron*) ;
  
  bool isBarrel_50(reco::GsfElectron*) ;
  bool isEndcap_50(reco::GsfElectron*) ;
  int detectorRegion_50(reco::GsfElectron*) ;
  
  enum DetectorRegion{ kNone , kBarrel , kEndcap , kGap , kForward } ; // kGap, kForward not yet used
};

class HEEPCut_41_Et                  : HEEPCutBase{
public:
  HEEPCut_41_Et(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_eta                 : HEEPCutBase{
public:
  HEEPCut_41_eta(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_EcalDriven          : HEEPCutBase{
public:
  HEEPCut_41_EcalDriven(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_dEtaIn              : HEEPCutBase{
public:
  HEEPCut_41_dEtaIn(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_dPhiIn              : HEEPCutBase{
public:
  HEEPCut_41_dPhiIn(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_HOverE              : HEEPCutBase{
public:
  HEEPCut_41_HOverE(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_SigmaIetaIeta       : HEEPCutBase{
public:
  HEEPCut_41_SigmaIetaIeta(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_E1x5OverE5x5        : HEEPCutBase{
public:
  HEEPCut_41_E1x5OverE5x5(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_E2x5OverE5x5        : HEEPCutBase{
public:
  HEEPCut_41_E2x5OverE5x5(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_isolEMHadDepth1     : HEEPCutBase{
private:
  float rho_ ;
  float EcalHcal1EffAreaBarrel_  ;
  float EcalHcal1EffAreaEndcaps_ ;
public:
  HEEPCut_41_isolEMHadDepth1(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setRho(float) ;
  void setEcalHcal1EffAreaBarrel (float) ;
  void setEcalHcal1EffAreaEndcaps(float) ;
} ;
class HEEPCut_41_IsolPtTrks          : HEEPCutBase{
public:
  HEEPCut_41_IsolPtTrks(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_missingHits         : HEEPCutBase{
public:
  HEEPCut_41_missingHits(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_41_dxyFirstPV          : HEEPCutBase{
private:
  math::XYZPoint firstPV_ ;
public:
  HEEPCut_41_dxyFirstPV(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setFirstPV(math::XYZPoint*) ;
} ;

class HEEPCut_50_50ns_Et             : HEEPCutBase{
public:
  HEEPCut_50_50ns_Et(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_eta            : HEEPCutBase{
public:
  HEEPCut_50_50ns_eta(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_EcalDriven     : HEEPCutBase{
public:
  HEEPCut_50_50ns_EcalDriven(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_dEtaIn         : HEEPCutBase{
public:
  HEEPCut_50_50ns_dEtaIn(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_dPhiIn         : HEEPCutBase{
public:
  HEEPCut_50_50ns_dPhiIn(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_HOverE         : HEEPCutBase{
public:
  HEEPCut_50_50ns_HOverE(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_SigmaIetaIeta  : HEEPCutBase{
public:
  HEEPCut_50_50ns_SigmaIetaIeta(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_E1x5OverE5x5   : HEEPCutBase{
public:
  HEEPCut_50_50ns_E1x5OverE5x5(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_E2x5OverE5x5   : HEEPCutBase{
public:
  HEEPCut_50_50ns_E2x5OverE5x5(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_isolEMHadDepth1: HEEPCutBase{
private:
  float rho_ ;
  float EcalHcal1EffAreaBarrel_  ;
  float EcalHcal1EffAreaEndcaps_ ;
public:
  HEEPCut_50_50ns_isolEMHadDepth1(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setRho(float) ;
  void setEcalHcal1EffAreaBarrel (float) ;
  void setEcalHcal1EffAreaEndcaps(float) ;
} ;
class HEEPCut_50_50ns_IsolPtTrks     : HEEPCutBase{
public:
  HEEPCut_50_50ns_IsolPtTrks(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_missingHits    : HEEPCutBase{
public:
  HEEPCut_50_50ns_missingHits(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_50ns_dxyFirstPV     : HEEPCutBase{
private:
  math::XYZPoint firstPV_ ;
public:
  HEEPCut_50_50ns_dxyFirstPV(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setFirstPV(math::XYZPoint*) ;
} ;

class HEEPCut_50_25ns_Et             : HEEPCutBase{
public:
  HEEPCut_50_25ns_Et(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_eta            : HEEPCutBase{
public:
  HEEPCut_50_25ns_eta(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_EcalDriven     : HEEPCutBase{
public:
  HEEPCut_50_25ns_EcalDriven(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_dEtaIn         : HEEPCutBase{
public:
  HEEPCut_50_25ns_dEtaIn(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_dPhiIn         : HEEPCutBase{
public:
  HEEPCut_50_25ns_dPhiIn(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_HOverE         : HEEPCutBase{
public:
  HEEPCut_50_25ns_HOverE(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_SigmaIetaIeta  : HEEPCutBase{
public:
  HEEPCut_50_25ns_SigmaIetaIeta(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_E1x5OverE5x5   : HEEPCutBase{
public:
  HEEPCut_50_25ns_E1x5OverE5x5(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_E2x5OverE5x5   : HEEPCutBase{
public:
  HEEPCut_50_25ns_E2x5OverE5x5(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_isolEMHadDepth1: HEEPCutBase{
private:
  float rho_ ;
  float EcalHcal1EffAreaBarrel_  ;
  float EcalHcal1EffAreaEndcaps_ ;
public:
  HEEPCut_50_25ns_isolEMHadDepth1(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setRho(float) ;
  void setEcalHcal1EffAreaBarrel (float) ;
  void setEcalHcal1EffAreaEndcaps(float) ;
} ;
class HEEPCut_50_25ns_IsolPtTrks     : HEEPCutBase{
public:
  HEEPCut_50_25ns_IsolPtTrks(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_missingHits    : HEEPCutBase{
public:
  HEEPCut_50_25ns_missingHits(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_50_25ns_dxyFirstPV     : HEEPCutBase{
private:
  math::XYZPoint firstPV_ ;
public:
  HEEPCut_50_25ns_dxyFirstPV(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setFirstPV(math::XYZPoint*) ;
} ;

class HEEPCutCollection{
private:
  bool status_ ;
  std::string name_ ;
  std::string branchName_n_ ;
  
  IIHEModuleHEEP* parent_mod_ ;
  
  std::vector<HEEPCutBase*      > listOfCuts_ ;
  std::vector<HEEPCutCollection*> listOfCutCollections_ ;
  
  // These variables keep track of the order of cuts/collections to preserve cutflow order
  std::vector<int> cutTypes_ ;
  int cutIndex_ ;
  int collectionIndex_ ;
  
  int nPass_ ;
  
  enum cutType{ kCut , kCollection } ;
public:
  HEEPCutCollection(std::string, IIHEModuleHEEP*) ;
  ~HEEPCutCollection() ;
  
  void beginEvent() ;
  void endEvent() ;
  
  void addCut(HEEPCutBase*) ;
  void addCutCollection(HEEPCutCollection*) ;
  bool applyCuts(reco::GsfElectron*, bool) ;
  bool applyCuts(reco::GsfElectron*) ;
  bool getStatus(){ return status_ ; }
  void makeCutflowHistogram() ;
  void config() ;
};

#endif
