#include "UserCode/IIHETree/interface/IIHEModuleMCTruth.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleMCTruth::IIHEModuleMCTruth(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
  pt_threshold_ = iConfig.getUntrackedParameter<double>("MCTruth_ptThreshold", 10.0) ;
  m_threshold_  = iConfig.getUntrackedParameter<double>("MCTruth_mThreshold" , 20.0) ;
}
IIHEModuleMCTruth::~IIHEModuleMCTruth(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleMCTruth::beginJob(){
  setBranchType(kVectorInt) ;
  addBranch("mc_index") ;
  addBranch("mc_pdgId") ;
  addBranch("mc_mother_index") ;
  addBranch("mc_mother_pdgId") ;
  addBranch("mc_charge") ;
  addBranch("mc_status") ;
  addBranch("mc_numberOfDaughters", kVectorUInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("mc_mass") ;
  addBranch("mc_pt") ;
  addBranch("mc_eta") ;
  addBranch("mc_phi") ;
  addBranch("mc_energy") ;
  
  whitelist_ = parent_->getMCTruthWhitelist() ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMCTruth::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Handle<GenParticleCollection> pGenParticles ;
  iEvent.getByLabel("genParticles", pGenParticles) ;
  GenParticleCollection genParticles(pGenParticles->begin(),pGenParticles->end());
  
  // These variables are used to match up mothers to daughters at the end
  std::vector<float> mother_pt  ;
  std::vector<float> mother_eta ;
  std::vector<float> mother_phi ;
  std::vector<float> mc_pt  ;
  std::vector<float> mc_eta ;
  std::vector<float> mc_phi ;
  int counter = 0 ;
  for(GenParticleCollection::const_iterator mc_iter=genParticles.begin() ; mc_iter!=genParticles.end() ; ++mc_iter){
    int pdgId = mc_iter->pdgId() ;
    float pt  = mc_iter->pt()  ;
    float eta = mc_iter->eta() ;
    float phi = mc_iter->phi() ;
    
    // First check the whitelist
    bool whitelist_accept = false ;
    for(unsigned int i=0 ; i<whitelist_.size() ; ++i){
      if(abs(pdgId)==abs(whitelist_.at(i))){
        whitelist_accept = true ;
        break ;
      }
    }
    
    // Ignore particles with exactly one daughter (X => X => X etc)
    bool daughters_accept = (mc_iter->numberOfDaughters()!=1) ;
    
    // Remove unphysical objects
    bool nonZeroPt_accept = (pt>1e-3) ;
    
    // Now check the thresholds
    bool thresholds_accept = false ;
    if(             pt>pt_threshold_) thresholds_accept = true  ;
    if(mc_iter->mass()> m_threshold_) thresholds_accept = true  ;
    
    // Now combine them all
    bool accept = (whitelist_accept && daughters_accept && thresholds_accept && nonZeroPt_accept) ;
    if(false==accept) continue ;
    
    const Candidate* mother = mc_iter->mother() ;
    while(abs(mother->pdgId())==abs(pdgId)){ mother = mother->mother() ; }
    
    mc_pt .push_back(pt ) ;
    mc_eta.push_back(eta) ;
    mc_phi.push_back(phi) ;
    mother_pt .push_back(mother->pt ()) ;
    mother_eta.push_back(mother->eta()) ;
    mother_phi.push_back(mother->phi()) ;
    
    store("mc_index" , counter          ) ;
    store("mc_pdgId" , pdgId            ) ;
    store("mc_charge", mc_iter->charge()) ;
    store("mc_status", mc_iter->status()) ;
    store("mc_mass"  , mc_iter->mass()  ) ;
    store("mc_pt"    , pt               ) ;
    store("mc_eta"   , eta              ) ;
    store("mc_phi"   , phi              ) ;
    store("mc_numberOfDaughters", (unsigned int)(mc_iter->numberOfDaughters())) ;
    store("mc_energy", mc_iter->energy()) ;
    store("mc_mother_pdgId", mother->pdgId()) ;
    counter++ ;
  }
  for(unsigned int i=0 ; i<mother_pt.size() ; ++i){
    int mother_index = -1 ;
    float best_DR = 1e6 ;
    for(unsigned int j=0 ; j<mc_pt.size() ; ++j){
      float DR = pow(mother_eta.at(i)-mc_eta.at(j), 2) + pow(mother_phi.at(i)-mc_phi.at(j), 2) ;
      if(DR<best_DR && DR<0.1){
        best_DR = DR ;
        mother_index = j ;
      }
    }
    store("mc_mother_index", mother_index) ;
  }
}

void IIHEModuleMCTruth::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMCTruth::beginEvent(){}
void IIHEModuleMCTruth::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMCTruth::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMCTruth);
