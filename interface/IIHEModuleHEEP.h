#ifndef UserCode_IIHETree_IIHEModuleHEEP_h
#define UserCode_IIHETree_IIHEModuleHEEP_h

class HEEPCut_41_isolEMHadDepth1 ;
class HEEPCut_41_dxyFirstPV ;
class HEEPCut_50_50ns_isolEMHadDepth1 ;
class HEEPCut_50_50ns_dxyFirstPV ;
class HEEPCut_50_25ns_isolEMHadDepth1 ;
class HEEPCut_50_25ns_dxyFirstPV ;
class HEEPCutCollection ;

#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "UserCode/IIHETree/interface/HEEPCut.h"

// class decleration
class IIHEModuleHEEP : public IIHEModule {
private:
  double EcalHcal1EffAreaBarrel_  ;
  double EcalHcal1EffAreaEndcaps_ ;
  double rho_ ;
  float triggerDeltaRThreshold_ ;
  
  HEEPCut_41_isolEMHadDepth1* cut_41_isolEMHadDepth1_ ;
  HEEPCut_41_dxyFirstPV*      cut_41_dxyFirstPV_ ;
  HEEPCutCollection* HEEPCutflow_41_ID_        ;
  HEEPCutCollection* HEEPCutflow_41_isolation_ ;
  HEEPCutCollection* HEEPCutflow_41_total_     ;
  
  HEEPCut_50_50ns_isolEMHadDepth1* cut_50_50ns_isolEMHadDepth1_ ;
  HEEPCut_50_50ns_dxyFirstPV*      cut_50_50ns_dxyFirstPV_ ;
  HEEPCutCollection* HEEPCutflow_50_50ns_ID_        ;
  HEEPCutCollection* HEEPCutflow_50_50ns_isolation_ ;
  HEEPCutCollection* HEEPCutflow_50_50ns_total_     ;
  
  HEEPCut_50_25ns_isolEMHadDepth1* cut_50_25ns_isolEMHadDepth1_ ;
  HEEPCut_50_25ns_dxyFirstPV*      cut_50_25ns_dxyFirstPV_ ;
  HEEPCutCollection* HEEPCutflow_50_25ns_ID_        ;
  HEEPCutCollection* HEEPCutflow_50_25ns_isolation_ ;
  HEEPCutCollection* HEEPCutflow_50_25ns_total_     ;
public:
  explicit IIHEModuleHEEP(const edm::ParameterSet& iConfig);
  ~IIHEModuleHEEP();
  
  void   pubBeginJob(){   beginJob() ; } ;
  void pubBeginEvent(){ beginEvent() ; } ;
  void   pubEndEvent(){   endEvent() ; } ;
  virtual void pubAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ analyze(iEvent, iSetup) ; } ;
  
  virtual void beginEvent() ;
  virtual void endEvent() ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  
  std::vector<HEEPCutCollection*> HEEPCutflows_ ;
};

#endif
