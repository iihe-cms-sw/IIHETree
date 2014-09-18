#ifndef UserCode_IIHETree_IIHEModuleTrigger_h
#define UserCode_IIHETree_IIHEModuleTrigger_h

#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "UserCode/IIHETree/interface/TriggerObject.h"

// class decleration
class IIHEModuleTrigger : public IIHEModule {
public:
  explicit IIHEModuleTrigger(const edm::ParameterSet& iConfig);
  ~IIHEModuleTrigger();
  
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
  
private:
  void addBranches() ;
  
  bool addL1TriggerElecton(std::string) ;
  bool addHLTriggerElecton(std::string, float) ;
  
  std::vector<L1Trigger*> L1TriggersElectron_ ;
  std::vector<HLTrigger*> HLTriggersElectron_ ;
  
  std::string branchPrefixElectronMatch_ ;
  
  std::vector<std::string>  L1FilterNamesPhoton_   ;
  std::vector<std::string>  L1FilterNamesElectron_ ;
  std::vector<std::string>  L1FilterNamesMuon_     ;
  
  std::vector<std::string> HLTFilterNamesPhoton_   ;
  std::vector<std::string> HLTFilterNamesElectron_ ;
  std::vector<std::string> HLTFilterNamesMuon_     ;
};
#endif
