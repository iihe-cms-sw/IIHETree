#ifndef UserCode_IIHETree_HEEPCut_h
#define UserCode_IIHETree_HEEPCut_h

#include <TH1.h>

#include "UserCode/IIHETree/interface/IIHEModuleHEEP.h"
// These classes are used to make it easier to keep track of HEEP cuts
// Each class corresponds to a different cut
// The cuts can be arranged into cut collections (eg HEEP41, HEEP50) to process a whole list of cuts at once

class IIHEModuleHEEP ;

class HEEPCutBase{
private:
  std::string name_ ;
  bool status_ ;
  int index_ ;
public:
  HEEPCutBase(std::string) ;
  ~HEEPCutBase() ;
  virtual bool applyCut(reco::GsfElectronCollection::const_iterator) ;
  bool addBranch(IIHEModuleHEEP*) ;
  void store(IIHEModuleHEEP*) ;
  void reset() ;
  void setStatus(bool) ;
  bool getStatus() ;
  std::string name(){ return name_ ; }
  
  bool isBarrel(reco::GsfElectronCollection::const_iterator) ;
  bool isEndcap(reco::GsfElectronCollection::const_iterator) ;
  int detectorRegion(reco::GsfElectronCollection::const_iterator) ;
  
  enum DetectorRegion{ kNone , kBarrel , kEndcap } ;
};

class HEEPCut_41_Et             : HEEPCutBase{
public:
  HEEPCut_41_Et(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_pt             : HEEPCutBase{
public:
  HEEPCut_41_pt(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_dEta           : HEEPCutBase{
public:
  HEEPCut_41_dEta(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_crack          : HEEPCutBase{
public:
  HEEPCut_41_crack(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_dEtaIn         : HEEPCutBase{
public:
  HEEPCut_41_dEtaIn(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_dPhiIn         : HEEPCutBase{
public:
  HEEPCut_41_dPhiIn(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_HOverE         : HEEPCutBase{
public:
  HEEPCut_41_HOverE(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_SigmaIetaIeta  : HEEPCutBase{
public:
  HEEPCut_41_SigmaIetaIeta(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_E2x5OverE5x5   : HEEPCutBase{
public:
  HEEPCut_41_E2x5OverE5x5(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_isolEMHadDepth1: HEEPCutBase{
private:
  float rho_ ;
  float EcalHcal1EffAreaBarrel_  ;
  float EcalHcal1EffAreaEndcaps_ ;
public:
  HEEPCut_41_isolEMHadDepth1(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
  void setRho(float) ;
  void setEcalHcal1EffAreaBarrel (float) ;
  void setEcalHcal1EffAreaEndcaps(float) ;
} ;
class HEEPCut_41_HcalIso2       : HEEPCutBase{
public:
  HEEPCut_41_HcalIso2(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_IsolPtTrks     : HEEPCutBase{
public:
  HEEPCut_41_IsolPtTrks(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_EcalDriven     : HEEPCutBase{
public:
  HEEPCut_41_EcalDriven(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_invalid        : HEEPCutBase{
public:
  HEEPCut_41_invalid(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_missingHits    : HEEPCutBase{
public:
  HEEPCut_41_missingHits(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_conversion     : HEEPCutBase{
public:
  HEEPCut_41_conversion(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
} ;
class HEEPCut_41_dxyFirstPV     : HEEPCutBase{
private:
  math::XYZPoint firstPV_ ;
public:
  HEEPCut_41_dxyFirstPV(std::string) ;
  bool applyCut(reco::GsfElectronCollection::const_iterator) ;
  void setFirstPV(math::XYZPoint) ;
} ;

class HEEPCutCollection{
private:
  bool status_ ;
  std::string name_ ;
  std::vector<HEEPCutBase*      > listOfCuts_ ;
  std::vector<HEEPCutCollection*> listOfCutCollections_ ;
  
  // These variables keep track of the order of cuts/collections to preserve cutflow order
  std::vector<int> cutTypes_ ;
  int cutIndex_ ;
  int collectionIndex_ ;
  
  // Histograms to keep track of things
  TH1F* hNElectrons_     ;
  TH1F* hNElectronsPass_ ;
  
  enum cutType{ kCut , kCollection } ;
public:
  HEEPCutCollection(std::string) ;
  ~HEEPCutCollection() ;
  void addCut(HEEPCutBase*) ;
  void addCutCollection(HEEPCutCollection*) ;
  bool applyCuts(reco::GsfElectronCollection::const_iterator, IIHEModuleHEEP*) ;
  bool getStatus(){ return status_ ; }
  void makeCutflowHistogram() ;
  void config(IIHEModuleHEEP*) ;
};

#endif
