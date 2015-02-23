#include "UserCode/IIHETree/interface/IIHEModuleSuperCluster.h"
#include "UserCode/IIHETree/interface/EtSort.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

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
  edm::Handle<reco::SuperClusterCollection> pHybridSuperClusters;
  edm::Handle<reco::SuperClusterCollection> pIslandSuperClusters;
  iEvent.getByLabel("correctedHybridSuperClusters"               ,"",pHybridSuperClusters);
  iEvent.getByLabel("correctedMulti5x5SuperClustersWithPreshower","",pIslandSuperClusters);
  const reco::SuperClusterCollection *hybridSuperClusters = pHybridSuperClusters.product() ;
  const reco::SuperClusterCollection *islandSuperClusters = pIslandSuperClusters.product() ;
  
  // Merge these two supercluster collections into one (sclusters collection)
  // There's probably a slicker way to do this...
  std::vector<const reco::SuperCluster*>    sclusters ;
  std::vector<reco::SuperClusterRef>     refsclusters ;
  for(reco::SuperClusterCollection::const_iterator hsc = hybridSuperClusters->begin() ; hsc!=hybridSuperClusters->end() ; hsc++ ){ sclusters.push_back(&(*hsc)) ; }
  for(reco::SuperClusterCollection::const_iterator isc = islandSuperClusters->begin() ; isc!=islandSuperClusters->end() ; isc++ ){ sclusters.push_back(&(*isc)) ; }
  for(unsigned int i=0 ; i<hybridSuperClusters->size() ; i++){ refsclusters.push_back(reco::SuperClusterRef(pHybridSuperClusters,i)) ; }
  for(unsigned int i=0 ; i<islandSuperClusters->size() ; i++){ refsclusters.push_back(reco::SuperClusterRef(pIslandSuperClusters,i)) ; }
  
  // Sort all the refSC and SC by transverse energy
  std::sort(refsclusters.begin(),refsclusters.end(),refScEtGreater) ;
  std::sort(   sclusters.begin(),   sclusters.end(),   scEtGreater) ;
  
  // Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();
  
  math::XYZPoint firstpvertex(0.0,0.0,0.0) ;
  float pv_z = -999 ;
  // We take only the first primary vertex, i.e. the one with the electrons
  if(pvcoll->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin() ;
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z()) ;
    pv_z = firstpv->z() ;
  }
  
  store("sc_n", (unsigned int) sclusters.size()) ;
  for(unsigned int i_sc=0 ; i_sc<sclusters.size() ; i_sc++){
    reco::SuperCluster* sc = (reco::SuperCluster*)sclusters.at(i_sc) ;
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
