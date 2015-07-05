#include "UserCode/IIHETree/interface/IIHEModuleLeptonsAccept.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleLeptonsAccept::IIHEModuleLeptonsAccept(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
}
IIHEModuleLeptonsAccept::~IIHEModuleLeptonsAccept(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleLeptonsAccept::beginJob(){
  nAcceptElEl_ = 0 ;
  nAcceptElMu_ = 0 ;
  nAcceptAll_  = 0 ;
}

// ------------ method called to for each event  ------------
void IIHEModuleLeptonsAccept::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  reco::MuonCollection muons = parent_->getMuonCollection() ;

  int nEl30 = 0 ;
  int nMu30 = 0 ;
  
  for(reco::GsfElectronCollection::const_iterator gsfiter=electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
    float pt = gsfiter->pt() ;
    float HEEP_ET  = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
    if(pt>30 || HEEP_ET>30) nEl30++ ;
  }
  for(reco::MuonCollection::const_iterator muiter = muons.begin(); muiter!=muons.end() ; ++muiter){
    float pt = muiter->pt() ;
    if(pt>30) nMu30++ ;
  }
  
  bool acceptElEl = (nEl30>=2) ;
  bool acceptElMu = (nEl30>=1 && nMu30>=1) ;
  bool acceptThisEvent = (acceptElEl || acceptElMu) ;
  
  // Save the event if we see something we like
  if(acceptThisEvent){
    acceptEvent() ;
    nAcceptAll_++ ;
  }
  if(acceptElEl) nAcceptElEl_ ++ ;
  if(acceptElMu) nAcceptElMu_ ++ ;
}

void IIHEModuleLeptonsAccept::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleLeptonsAccept::beginEvent(){}
void IIHEModuleLeptonsAccept::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleLeptonsAccept::endJob(){
  std::cout << std::endl << "IIHEModuleLeptonsAccept report:" << std::endl ;
  std::cout << "  nAcceptElEl  = " << nAcceptElEl_ << std::endl ;
  std::cout << "  nAcceptMuMu  = " << nAcceptElMu_ << std::endl ;
  std::cout << "  nAcceptAll  = " << nAcceptAll_   << std::endl ;
}

DEFINE_FWK_MODULE(IIHEModuleLeptonsAccept);
