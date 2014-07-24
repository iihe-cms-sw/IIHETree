#include "UserCode/IIHETree/interface/IIHEModuleVertex.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;

IIHEModuleVertex::IIHEModuleVertex(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
}
IIHEModuleVertex::~IIHEModuleVertex(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleVertex::beginJob(){
  addBranch("pv_n", kInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("pv_x") ;
  addBranch("pv_y") ;
  addBranch("pv_z") ;
  addBranch("pv_isValid", kVectorBool) ;
  addBranch("pv_normalizedChi2", kVectorFloat) ;
  setBranchType(kVectorInt) ;
  addBranch("pv_ndof") ;
  addBranch("pv_nTracks") ;
  addBranch("pv_totTrackSize") ;
}

// ------------ method called to for each event  ------------
void IIHEModuleVertex::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // Get the beamspot from the Event:
  // The beamspot is passed to the IIHEAnalysis class, so we call it from parent_
  // Don't forget to declare IIHEModuleVertex as a friend of IIHEAnalysis!
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(parent_->beamSpotLabel_, theBeamSpot);

  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint firstpvertex      (0.0,0.0,0.0) ;
  math::XYZPoint firstpvertexwithBS(0.0,0.0,0.0) ;

  // Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();
  
  Handle<reco::VertexCollection> primaryVertexCollwithBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVertexCollwithBS);
  const reco::VertexCollection* pvcollwithBS = primaryVertexCollwithBS.product();
  
  if(pvcoll->size()>0){
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  if(pvcollwithBS->size()>0){
    reco::VertexCollection::const_iterator firstpv = pvcollwithBS->begin();
    firstpvertexwithBS.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }

  store("pv_n", (int)(pvcoll->size())) ;
  for(reco::VertexCollection::const_iterator pvIt = pvcoll->begin(); pvIt != pvcoll->end(); ++pvIt){
    store("pv_x"             , pvIt->x()) ;
    store("pv_y"             , pvIt->y()) ;   
    store("pv_z"             , pvIt->z()) ;  
    store("pv_isValid"       , pvIt->isValid()) ;
    store("pv_ndof"          , (int)pvIt->ndof()) ;
    store("pv_nTracks"       , (int)(pvIt->nTracks())) ;
    store("pv_normalizedChi2", pvIt->normalizedChi2()) ;
    store("pv_totTrackSize"  , (int)(pvIt->tracksSize())) ;
  }
}

void IIHEModuleVertex::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleVertex::beginEvent(){}
void IIHEModuleVertex::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleVertex::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleVertex);
