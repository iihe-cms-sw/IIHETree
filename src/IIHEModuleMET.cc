#include "UserCode/IIHETree/interface/IIHEModuleMET.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h" 
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleMET::IIHEModuleMET(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleMET::~IIHEModuleMET(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleMET::beginJob(){
  setBranchType(kFloat) ;
  addBranch("MET_caloMet_et" ) ;
  addBranch("MET_caloMet_phi") ;
  addBranch("MET_pfMet_et"   ) ;
  addBranch("MET_pfMet_phi"  ) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMET::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<CaloMETCollection> pCaloMET;
  bool caloMetIsValid = iEvent.getByLabel("caloMet", pCaloMET);
  const CaloMETCollection *caloMET = pCaloMET.product();

  edm::Handle<PFMETCollection> pPFMET;
  bool pfMetIsValid = iEvent.getByLabel("pfMet", pPFMET);
  const PFMETCollection *PFMET = pPFMET.product();
  
  if(caloMetIsValid){
    store("MET_caloMet_et" , caloMET->begin()->et()  ) ;
    store("MET_caloMet_phi", caloMET->begin()->phi() ) ;
  }
  if(pfMetIsValid){
    store("MET_pfMet_et"   , PFMET->begin()->et()    ) ;
    store("MET_pfMet_phi"  , PFMET->begin()->phi()   ) ;
  }
}

void IIHEModuleMET::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMET::beginEvent(){}
void IIHEModuleMET::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMET::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMET);
