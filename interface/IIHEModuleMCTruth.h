#ifndef UserCode_IIHETree_IIHEModuleMCTruth_h
#define UserCode_IIHETree_IIHEModuleMCTruth_h

#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "UserCode/IIHETree/interface/MCTruthObject.h"

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"

// class decleration
class IIHEModuleMCTruth : public IIHEModule {
public:
  explicit IIHEModuleMCTruth(const edm::ParameterSet& iConfig);
  ~IIHEModuleMCTruth();
  
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
  std::vector<int> whitelist_ ;
  double pt_threshold_ ;
  double  m_threshold_ ;
};
#endif
