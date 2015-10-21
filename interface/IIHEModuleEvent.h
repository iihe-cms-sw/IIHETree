#ifndef UserCode_IIHETree_IIHEModuleEvent_h
#define UserCode_IIHETree_IIHEModuleEvent_h

#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "DataFormats/JetReco/interface/Jet.h"

// class decleration
class IIHEModuleEvent : public IIHEModule {
private:
  edm::InputTag rhoLabel_ ;

public:
  explicit IIHEModuleEvent(const edm::ParameterSet& iConfig);
  ~IIHEModuleEvent();
  
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
};
#endif
