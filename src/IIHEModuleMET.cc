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
  addBranch("MET_met_et") ;
  addBranch("MET_met_phi") ;
  
  //addBranch("MET_htMetKT4_et") ;
  //addBranch("MET_htMetKT4_phi") ;
  
  addBranch("MET_pfMet_et") ;
  addBranch("MET_pfMet_phi") ;
  
  addBranch("pfType1CorrectedMet_met_et") ;
  addBranch("pfType1CorrectedMet_met_phi") ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMET::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<CaloMETCollection> pCaloMET;
  bool calometisvalid = iEvent.getByLabel("met", pCaloMET);
  const CaloMETCollection *caloMET  = pCaloMET.product();

  //edm::Handle<METCollection> pMET;
  //bool metisvalid = iEvent.getByLabel("htMetKT4", pMET);
  //const METCollection *MET  = pMET.product();

  edm::Handle<PFMETCollection> pPFMET;
  bool pfmetisvalid = iEvent.getByLabel("pfMet", pPFMET);
  const PFMETCollection *PFMET  = pPFMET.product();

  edm::Handle<PFMETCollection> pPFMETcor;
  bool pfmetcorisvalid = iEvent.getByLabel("pfType1CorrectedMet", pPFMETcor);
  const PFMETCollection *PFMETcor  = pPFMETcor.product();
  
  if(calometisvalid){
    store("MET_met_et"                 , caloMET->begin()->et()  ) ;
    store("MET_met_phi"                , caloMET->begin()->phi() ) ;
  }
  //if(metisvalid){
  //  store("MET_htMetKT4_et"            , MET->begin()->et()      ) ;
  //  store("MET_htMetKT4_phi"           , MET->begin()->phi()     ) ;
  //}
  if(pfmetisvalid){
    store("MET_pfMet_et"               , PFMET->begin()->et()    ) ;
    store("MET_pfMet_phi"              , PFMET->begin()->phi()   ) ;
  }
  if(pfmetcorisvalid){
    store("pfType1CorrectedMet_met_et" , PFMETcor->begin()->et() ) ;
    store("pfType1CorrectedMet_met_phi", PFMETcor->begin()->phi()) ;
  }
}

void IIHEModuleMET::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMET::beginEvent(){}
void IIHEModuleMET::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMET::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMET);
