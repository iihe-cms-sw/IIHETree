#include "UserCode/IIHETree/interface/IIHEModuleSuperCluster.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleSuperCluster::IIHEModuleSuperCluster(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleSuperCluster::~IIHEModuleSuperCluster(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleSuperCluster::beginJob(){
  addBranch("sc_n",kUInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("sc_energy") ;
  addBranch("sc_eta") ;
  addBranch("sc_etacorr") ;
  addBranch("sc_theta") ;
  addBranch("sc_thetacorr") ;
  addBranch("sc_et") ;
  addBranch("sc_phi") ;
  addBranch("sc_px") ;
  addBranch("sc_py") ;
  addBranch("sc_pz") ;
  addBranch("sc_x") ;
  addBranch("sc_y") ;
  addBranch("sc_z") ;
}

// ------------ method called to for each event  ------------
void IIHEModuleSuperCluster::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  std::vector<const reco::SuperCluster*> superclusters = parent_->getSuperClusters() ;
  
  // Retrieve primary vertex collection
  math::XYZPoint* firstPV = parent_->getFirstPrimaryVertex() ;
  float pv_z = firstPV->z() ;
  
  store("sc_n", (unsigned int) superclusters.size()) ;
  for(unsigned int i_sc=0 ; i_sc<superclusters.size() ; i_sc++){
    reco::SuperCluster* sc = (reco::SuperCluster*) superclusters.at(i_sc) ;
    float sc_energy = sc->rawEnergy()+sc->preshowerEnergy() ;
    float sc_et     = sc_energy/cosh(sc->eta()) ;
      
    store("sc_eta"      , sc->eta()) ;
    store("sc_etacorr"  , etacorr( sc->eta(), pv_z, sc->position().z() )) ;
    store("sc_theta"    , 2.*atan(exp(-1.*sc->eta()))) ;
    store("sc_thetacorr", 2.*atan(exp(-1.*etacorr( sc->eta(), pv_z, sc->position().z() ) ))) ;
    store("sc_phi"      , sc->phi()) ;
    store("sc_energy"   , sc_energy) ;
    store("sc_et"       , sc_et    );
    store("sc_px"       , sc_et*cos(sc->phi())) ;
    store("sc_py"       , sc_et*sin(sc->phi()));
    store("sc_pz"       , sc_energy*tanh(sc->eta())) ;
    store("sc_x"        , sc->position().x()) ;
    store("sc_y"        , sc->position().y()) ;
    store("sc_z"        , sc->position().z()) ;
  }
}

void IIHEModuleSuperCluster::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleSuperCluster::beginEvent(){}
void IIHEModuleSuperCluster::endEvent(){}

// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleSuperCluster::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleSuperCluster);
