#include "UserCode/IIHETree/interface/IIHEModuleEvent.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleEvent::IIHEModuleEvent(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleEvent::~IIHEModuleEvent(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleEvent::beginJob(){
  setBranchType(kUInt) ;
  addBranch("ev_event"          ) ;
  addBranch("ev_run"            ) ;
  addBranch("ev_luminosityBlock") ;
  
  addBranch("ev_rho", kFloat) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleEvent::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  store("ev_event"          , ((unsigned int) (iEvent.id().event()          ))) ;
  store("ev_run"            , ((unsigned int) (iEvent.id().run()            ))) ;
  store("ev_luminosityBlock", ((unsigned int) (iEvent.id().luminosityBlock()))) ;
  
  edm::Handle<double> rhoHandle ;
  iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle) ;
  float rho = *rhoHandle ;
  store("ev_rho", rho) ;
}

void IIHEModuleEvent::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleEvent::beginEvent(){}
void IIHEModuleEvent::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleEvent::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleEvent);
