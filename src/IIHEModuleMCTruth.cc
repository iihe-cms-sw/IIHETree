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
  addBranch("mc_n", kUInt) ;
  setBranchType(kVectorInt) ;
  addBranch("mc_index") ;
  addBranch("mc_pdgId") ;
  addBranch("mc_charge") ;
  addBranch("mc_status") ;
  setBranchType(kVectorFloat) ;
  addBranch("mc_mass") ;
  addBranch("mc_px") ;
  addBranch("mc_py") ;
  addBranch("mc_pz") ;
  addBranch("mc_pt") ;
  addBranch("mc_eta") ;
  addBranch("mc_phi") ;
  addBranch("mc_energy") ;
  setBranchType(kVectorUInt) ;
  addBranch("mc_numberOfDaughters") ;
  addBranch("mc_numberOfMothers"  ) ;
  setBranchType(kVectorVectorInt) ;
  addBranch("mc_mother_index") ;
  addBranch("mc_mother_pdgId") ;
  setBranchType(kVectorVectorFloat) ;
  addBranch("mc_mother_px"    ) ;
  addBranch("mc_mother_py"    ) ;
  addBranch("mc_mother_pz"    ) ;
  addBranch("mc_mother_pt"    ) ;
  addBranch("mc_mother_eta"   ) ;
  addBranch("mc_mother_phi"   ) ;
  addBranch("mc_mother_energy") ;
  addBranch("mc_mother_mass"  ) ;
  
  whitelist_ = parent_->getMCTruthWhitelist() ;
}

// ------------ method called to for each event  ------------
void IIHEModuleMCTruth::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  Handle<GenParticleCollection> pGenParticles ;
  iEvent.getByLabel("genParticles", pGenParticles) ;
  GenParticleCollection genParticles(pGenParticles->begin(),pGenParticles->end());
  
  // These variables are used to match up mothers to daughters at the end
  int counter = 0 ;
  
  std::vector<MCTruthObject*> MCTruthRecord ;
  for(GenParticleCollection::const_iterator mc_iter=genParticles.begin() ; mc_iter!=genParticles.end() ; ++mc_iter){
    int pdgId = mc_iter->pdgId() ;
    float pt  = mc_iter->pt()  ;
    
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
    
    if(mc_iter->numberOfDaughters()==2 && abs(mc_iter->pdgId())==11 && false){
      // Just debugging to see how often a photon gets radiated
      const reco::Candidate* d1 = mc_iter->daughter(0) ;
      const reco::Candidate* d2 = mc_iter->daughter(1) ;
      std::cout << mc_iter->pdgId() << " " << mc_iter->px() << " " << mc_iter->py() << " " << mc_iter->pz() << " " << mc_iter->energy() << " " << mc_iter->mass() << std::endl ;
      std::cout << d1->pdgId() << " " << d1->px() << " " << d1->py() << " " << d1->pz() << " " << d1->energy() << " " << d1->mass() << std::endl ;
      std::cout << d2->pdgId() << " " << d2->px() << " " << d2->py() << " " << d2->pz() << " " << d2->energy() << " " << d2->mass() << std::endl ;
      std::cout << std::endl ;
    }
    
    // Remove unphysical objects
    bool nonZeroPt_accept = (pt>1e-3) ;
    
    // Now check the thresholds
    bool thresholds_accept = false ;
    if(             pt>pt_threshold_) thresholds_accept = true  ;
    if(mc_iter->mass()> m_threshold_) thresholds_accept = true  ;
    
    // Now combine them all
    bool accept = (whitelist_accept && daughters_accept && thresholds_accept && nonZeroPt_accept) ;
    if(false==accept) continue ;
    
    // Now go up the ancestry until we find the real parent
    const Candidate* parent = mc_iter->mother() ;
    const Candidate* child  = mc_iter->clone()  ;
    while(parent->pdgId()==pdgId){
      child  = parent ;
      parent = parent->mother() ;
    }
    
    // Create a truth record instance
    MCTruthObject* MCTruth = new MCTruthObject((reco::Candidate*)&*mc_iter) ;
    
    // Add all the mothers
    for(unsigned int mother_iter=0 ; mother_iter<child->numberOfMothers() ; ++mother_iter){
      MCTruth->addMother(child->mother(mother_iter)) ;
    }
    
    // Then push back the MC truth information
    MCTruthRecord.push_back(MCTruth) ;
    counter++ ;
  }
  for(unsigned int i=0 ; i<MCTruthRecord.size() ; ++i){
    MCTruthObject* ob = MCTruthRecord.at(i) ;
    std::vector<int  > mc_mother_index ;
    std::vector<int  > mc_mother_pdgId ;
    std::vector<float> mc_mother_px ;
    std::vector<float> mc_mother_py ;
    std::vector<float> mc_mother_pz ;
    std::vector<float> mc_mother_pt ;
    std::vector<float> mc_mother_eta ;
    std::vector<float> mc_mother_phi ;
    std::vector<float> mc_mother_energy ;
    std::vector<float> mc_mother_mass ;
    for(unsigned int j=0 ; j<ob->nMothers() ; ++j){
      const reco::Candidate* mother = ob->getMother(j) ;
      if(mother){
        int mother_index_tmp = ob->matchMother(MCTruthRecord, j) ;
        mc_mother_index .push_back(mother_index_tmp) ;
        mc_mother_pdgId .push_back(mother->pdgId() ) ;
        mc_mother_px    .push_back(mother->px()    ) ;
        mc_mother_py    .push_back(mother->py()    ) ;
        mc_mother_pz    .push_back(mother->pz()    ) ;
        mc_mother_pt    .push_back(mother->pt()    ) ;
        mc_mother_eta   .push_back(mother->eta()   ) ;
        mc_mother_phi   .push_back(mother->phi()   ) ;
        mc_mother_energy.push_back(mother->energy()) ;
        mc_mother_mass  .push_back(mother->mass()  ) ;
      }
    }
    if(mc_mother_index.size()==0) mc_mother_index.push_back(0) ;
    store("mc_mother_index" , mc_mother_index ) ;
    store("mc_mother_pdgId" , mc_mother_pdgId ) ;
    store("mc_mother_px"    , mc_mother_px    ) ;
    store("mc_mother_py"    , mc_mother_py    ) ;
    store("mc_mother_pz"    , mc_mother_pz    ) ;
    store("mc_mother_pt"    , mc_mother_pt    ) ;
    store("mc_mother_eta"   , mc_mother_eta   ) ;
    store("mc_mother_phi"   , mc_mother_phi   ) ;
    store("mc_mother_energy", mc_mother_energy) ;
    store("mc_mother_mass"  , mc_mother_mass  ) ;
    
    store("mc_numberOfDaughters", (unsigned int)(ob->getCandidate()->numberOfDaughters())) ;
    store("mc_numberOfMothers"  , (unsigned int)(ob->nMothers())) ;
    
    store("mc_px"     , ob->getCandidate()->px()    ) ;
    store("mc_py"     , ob->getCandidate()->py()    ) ;
    store("mc_pz"     , ob->getCandidate()->pz()    ) ;
    store("mc_pt"     , ob->getCandidate()->pt()    ) ;
    store("mc_eta"    , ob->getCandidate()->eta()   ) ;
    store("mc_phi"    , ob->getCandidate()->phi()   ) ;
    store("mc_energy" , ob->getCandidate()->energy()) ;
    store("mc_mass"   , ob->getCandidate()->mass()  ) ;
    
    store("mc_index"  , i ) ;
    store("mc_pdgId"  , ob->getCandidate()->pdgId() ) ;
    store("mc_charge" , ob->getCandidate()->charge()) ;
    store("mc_status" , ob->getCandidate()->status()) ;
  }
  store("mc_n", (unsigned int)(MCTruthRecord.size())) ;
}

void IIHEModuleMCTruth::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleMCTruth::beginEvent(){}
void IIHEModuleMCTruth::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleMCTruth::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleMCTruth);
