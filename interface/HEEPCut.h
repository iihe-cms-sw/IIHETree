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
  float value_ ; // Value of the variable
  bool status_ ;
  int index_ ;
  int nPass_ ;
  int nPassCumulative_ ;
  IIHEModuleHEEP* parent_mod_ ;
  
  float barrelEtaUpper_ ;
  float endcapEtaLower_ ;
  float endcapEtaUpper_ ;
  
public:
  HEEPCutBase(std::string, IIHEModuleHEEP*) ;
  ~HEEPCutBase() ;
  virtual bool applyCut(reco::GsfElectron*, bool) ;
  bool addBranches() ;
  void store() ;
  void reset() ;
  void setStatus(bool, bool) ;
  bool getStatus() ;
  void setValue(float value){ value_ = value ; }
  float value(){ return value_ ; }
  std::string name(){ return name_ ; }
  
  void beginEvent() ;
  void endEvent() ;
  
  bool isBarrel(reco::GsfElectron*) ;
  bool isEndcap(reco::GsfElectron*) ;
  int detectorRegion(reco::GsfElectron*) ;
  
  void setBarrelLimits(float barrelEtaUpper){ barrelEtaUpper_ = barrelEtaUpper ; }
  void setEndcapLimits(float endcapEtaLower, float endcapEtaUpper){
    endcapEtaLower_ = endcapEtaLower ;
    endcapEtaUpper_ = endcapEtaUpper ;
  }
  
  enum DetectorRegion{ kNone , kBarrel , kEndcap , kGap , kForward } ; // kGap, kForward not yet used
};

class HEEPCut_Et: HEEPCutBase{
public:
  HEEPCut_Et(std::string, IIHEModuleHEEP*, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float thresholdBarrel_ ;
  float thresholdEndcap_ ;
} ;
class HEEPCut_eta: HEEPCutBase{
public:
  HEEPCut_eta(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_EcalDriven: HEEPCutBase{
public:
  HEEPCut_EcalDriven(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_dPhiIn: HEEPCutBase{
public:
  HEEPCut_dPhiIn(std::string, IIHEModuleHEEP*, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float thresholdBarrel_ ;
  float thresholdEndcap_ ;
} ;
class HEEPCut_SigmaIetaIeta: HEEPCutBase{
public:
  HEEPCut_SigmaIetaIeta(std::string, IIHEModuleHEEP*, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float threshold_ ;
} ;
class HEEPCut_E1x5OverE5x5: HEEPCutBase{
public:
  HEEPCut_E1x5OverE5x5(std::string, IIHEModuleHEEP*) ;
  bool applyCut(reco::GsfElectron*, bool) ;
} ;
class HEEPCut_E2x5OverE5x5: HEEPCutBase{
public:
  HEEPCut_E2x5OverE5x5(std::string, IIHEModuleHEEP*, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float thresholdE1x5_ ;
  float thresholdE2x5_ ;
} ;
class HEEPCut_isolEMHadDepth1: HEEPCutBase{
private:
  float rho_ ;
  float EcalHcal1EffAreaBarrel_  ;
  float EcalHcal1EffAreaEndcaps_ ;
  float constantTermBarrel_ ;
  float constantTermEndcapLowEt_ ;
  float constantTermEndcapHighEt_ ;
  float linearTermBarrel_ ;
  float linearTermEndcap_ ;
  float offsetTermEndcap_ ;
public:
  HEEPCut_isolEMHadDepth1(std::string, IIHEModuleHEEP*, float, float, float, float, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setRho(float) ;
  void setEcalHcal1EffAreaBarrel (float) ;
  void setEcalHcal1EffAreaEndcaps(float) ;
private:
} ;
class HEEPCut_IsolPtTrks: HEEPCutBase{
public:
  HEEPCut_IsolPtTrks(std::string, IIHEModuleHEEP*, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float thresholdBarrel_ ;
  float thresholdEndcap_ ;
} ;
class HEEPCut_missingHits: HEEPCutBase{
public:
  HEEPCut_missingHits(std::string, IIHEModuleHEEP*, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float threshold_ ;
} ;
class HEEPCut_dxyFirstPV: HEEPCutBase{
private:
  math::XYZPoint firstPV_ ;
  float thresholdBarrel_ ;
  float thresholdEndcap_ ;
public:
  HEEPCut_dxyFirstPV(std::string, IIHEModuleHEEP*, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
  void setFirstPV(math::XYZPoint*) ;
} ;


class HEEPCut_41_dEtaIn: HEEPCutBase{
public:
  HEEPCut_41_dEtaIn(std::string, IIHEModuleHEEP*, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float thresholdBarrel_ ;
  float thresholdEndcap_ ;
} ;
class HEEPCut_41_HOverE: HEEPCutBase{
public:
  HEEPCut_41_HOverE(std::string, IIHEModuleHEEP*, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float thresholdBarrel_ ;
  float thresholdEndcap_ ;
} ;


class HEEPCut_50_50ns_dEtaIn: HEEPCutBase{
public:
  HEEPCut_50_50ns_dEtaIn(std::string, IIHEModuleHEEP*, float, float, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float constantTermBarrel_ ;
  float linearTermBarrel_   ;
  float cutoffTermBarrel_   ;
  float thresholdEndcap_    ;
} ;
class HEEPCut_50_50ns_HOverE: HEEPCutBase{
public:
  HEEPCut_50_50ns_HOverE(std::string, IIHEModuleHEEP*, float, float, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float reciprocalTermBarrel_ ;
  float reciprocalTermEndcap_ ;
  float constantTermBarrel_ ;
  float constantTermEndcap_ ;
} ;


class HEEPCut_50_25ns_dEtaIn: HEEPCutBase{
public:
  HEEPCut_50_25ns_dEtaIn(std::string, IIHEModuleHEEP*, float, float, float, float, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float constantTermBarrel_ ;
  float linearTermBarrel_   ;
  float cutoffTermBarrel_   ;
  float constantTermEndcap_ ;
  float linearTermEndcap_   ;
  float cutoffTermEndcap_   ;
} ;
class HEEPCut_50_25ns_HOverE: HEEPCutBase{
public:
  HEEPCut_50_25ns_HOverE(std::string, IIHEModuleHEEP*, float, float, float, float) ;
  bool applyCut(reco::GsfElectron*, bool) ;
private:
  float reciprocalTermBarrel_ ;
  float reciprocalTermEndcap_ ;
  float constantTermBarrel_ ;
  float constantTermEndcap_ ;
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
  void config(float, float, float) ;
};

#endif
