#include "UserCode/IIHETree/interface/IIHEModuleHEEP.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"


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

IIHEModuleHEEP::IIHEModuleHEEP(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
  EcalHcal1EffAreaBarrel_  = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaBarrel" , 0.) ;
  EcalHcal1EffAreaEndcaps_ = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaEndcaps", 0.) ;
  rho_ = iConfig.getUntrackedParameter<double>("kt6PFJets:rho", 0.) ;
}
IIHEModuleHEEP::~IIHEModuleHEEP(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleHEEP::beginJob(){
  std::vector<int> MCPdgIdsToSave ;
  MCPdgIdsToSave.push_back(11) ; // Electron
  MCPdgIdsToSave.push_back(13) ; // Muon
  MCPdgIdsToSave.push_back(15) ; // Tau
  MCPdgIdsToSave.push_back(22) ; // Photon
  MCPdgIdsToSave.push_back(23) ; // Z boson
  MCPdgIdsToSave.push_back(24) ; // W boson
  MCPdgIdsToSave.push_back(25) ; // BEH boson
  MCPdgIdsToSave.push_back( 5) ; // b quark
  MCPdgIdsToSave.push_back( 6) ; // t quark
  MCPdgIdsToSave.push_back(32) ; // Z'  boson
  MCPdgIdsToSave.push_back(33) ; // Z'' boson
  MCPdgIdsToSave.push_back(34) ; // W'  boson
  addToMCTruthWhitelist(MCPdgIdsToSave) ;
  
  addBranch("HEEP_nHEEP", kUInt);
  
  // Preshower information
  setBranchType(kVectorFloat) ;
  addBranch("HEEP_eseffsixix") ;
  addBranch("HEEP_eseffsiyiy") ;
  addBranch("HEEP_eseffsirir") ;
  addBranch("HEEP_preshowerEnergy") ;
  
  addBranch("HEEP_e1x3") ;
  
  // Crystal information
  setBranchType(kVectorVectorFloat) ;
  addBranch("HEEP_crystal_energy") ;
  addBranch("HEEP_crystal_eta"   ) ;
  addBranch("HEEP_eshitsixix") ;
  addBranch("HEEP_eshitsiyiy") ;
  setBranchType(kVectorVectorInt) ;
  addBranch("HEEP_crystal_ietaorix") ;
  addBranch("HEEP_crystal_iphioriy") ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                     HEEP cutflow                                   //
  ////////////////////////////////////////////////////////////////////////////////////////
    
  HEEPCutflow_41_ID_        = new HEEPCutCollection("HEEP_cutflow41_ID"       ) ;
  HEEPCutflow_41_isolation_ = new HEEPCutCollection("HEEP_cutflow41_isolation") ;
  HEEPCutflow_41_total_     = new HEEPCutCollection("HEEP_cutflow41_total"    ) ;
  
  // These cuts must be updated so declare them separately
  cut_41_isolEMHadDepth1_ = new HEEPCut_41_isolEMHadDepth1("HEEPCut_cutflow41_isolEMHadDepth1") ;
  cut_41_dxyFirstPV_      = new HEEPCut_41_dxyFirstPV     ("HEEP_cutflow41_dxyFirstPV "       ) ;
  
  cut_41_isolEMHadDepth1_->setEcalHcal1EffAreaBarrel (EcalHcal1EffAreaBarrel_ ) ;
  cut_41_isolEMHadDepth1_->setEcalHcal1EffAreaEndcaps(EcalHcal1EffAreaEndcaps_) ;
  cut_41_isolEMHadDepth1_->HEEPCut_41_isolEMHadDepth1::setRho(rho_) ;
  
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_dEtaIn       ("HEEP_cutflow41_dEtaIn"       )) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_dPhiIn       ("HEEP_cutflow41_dPhiIn"       )) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_HOverE       ("HEEP_cutflow41_HOverE"       )) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_SigmaIetaIeta("HEEP_cutflow41_SigmaIetaIeta")) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_E2x5OverE5x5 ("HEEP_cutflow41_E2x5OverE5x5" )) ;
  
  HEEPCutflow_41_isolation_->addCut( (HEEPCutBase*) cut_41_isolEMHadDepth1_) ;
  HEEPCutflow_41_isolation_->addCut( (HEEPCutBase*) new HEEPCut_41_HcalIso2  ("HEEP_cutflow41_HcalIso2"  )) ;
  HEEPCutflow_41_isolation_->addCut( (HEEPCutBase*) new HEEPCut_41_IsolPtTrks("HEEP_cutflow41_IsolPtTrks")) ;
  
  HEEPCutflow_41_total_->addCutCollection(HEEPCutflow_41_ID_       ) ;
  HEEPCutflow_41_total_->addCutCollection(HEEPCutflow_41_isolation_) ;
  
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_Et         ("HEEP_cutflow41_Et"         )) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_pt         ("HEEP_cutflow41_pt"         )) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_dEta       ("HEEP_cutflow41_dEta"       )) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_crack      ("HEEP_cutflow41_crack"      )) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_EcalDriven ("HEEP_cutflow41_EcalDriven" )) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_invalid    ("HEEP_cutflow41_invalid"    )) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_missingHits("HEEP_cutflow41_missingHits")) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) new HEEPCut_41_conversion ("HEEP_cutflow41_conversion" )) ;
  HEEPCutflow_41_total_->addCut( (HEEPCutBase*) cut_41_dxyFirstPV_) ;
  
  HEEPCutflow_41_ID_       ->config(this) ;
  HEEPCutflow_41_isolation_->config(this) ;
  HEEPCutflow_41_total_    ->config(this) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleHEEP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
  iEvent.getByLabel("gedGsfElectrons","",pGsfElectrons);
  reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());
  
  // Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();
  
  math::XYZPoint firstpvertex(0.0,0.0,0.0) ;
  // We take only the first primary vertex, i.e. the one with the electrons
  if(pvcoll->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  cut_41_dxyFirstPV_->setFirstPV(firstpvertex) ;
  
  // Information for preshower
  edm::ESHandle<CaloGeometry> pGeometry ;
  iSetup.get<CaloGeometryRecord>().get(pGeometry) ;
  CaloGeometry* geometry = (CaloGeometry*) pGeometry.product() ;
  const CaloSubdetectorGeometry* geometryES = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower) ;
  CaloSubdetectorTopology* topology_p = 0 ;
  if(geometryES) topology_p = new EcalPreshowerTopology(geometry) ;
  
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  
  Handle<EcalRecHitCollection> EBhits;
  Handle<EcalRecHitCollection> EEhits;
  iEvent.getByLabel("reducedEcalRecHitsEB",EBhits);
  iEvent.getByLabel("reducedEcalRecHitsEE",EEhits);
  EcalClusterLazyTools lazytool(iEvent,iSetup,InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"),InputTag("reducedEcalRecHitsES"));
  
  unsigned int nHeepEle = 0 ;
  for(reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin() ; gsfiter!=gsfelectrons.end() ; ++gsfiter){
    // Required for preshower variables
    reco::SuperClusterRef cl_ref = gsfiter->superCluster();
    const reco::CaloClusterPtr seed = gsfiter->superCluster()->seed();
    
    HEEPCutflow_41_total_->applyCuts(gsfiter, this) ;
    
    // Preshower information
    // Get the preshower hits
    double x = gsfiter->superCluster()->x() ;
    double y = gsfiter->superCluster()->y() ;
    double z = gsfiter->superCluster()->z() ;
    store("HEEP_eshitsixix", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 1)) ;
    store("HEEP_eshitsiyiy", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 2)) ;
    store("HEEP_preshowerEnergy", gsfiter->superCluster()->preshowerEnergy()) ;
    
    store("HEEP_eseffsixix", lazytool.eseffsixix(*cl_ref)) ;
    store("HEEP_eseffsiyiy", lazytool.eseffsiyiy(*cl_ref)) ;
    store("HEEP_eseffsirir", lazytool.eseffsirir(*cl_ref)) ;
    store("HEEP_e1x3"      , lazytool.e1x3(*seed)        ) ;
    
    // Try to add info about rechit in the SC 
    // Strongly inspired from : http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/DaveC/src/printPhoton.cc
    //Crystal variables
    std::vector<float> gsf_crystal_energy  ;
    std::vector<int  > gsf_crystal_ietaorix;
    std::vector<int  > gsf_crystal_iphioriy;
    std::vector<float> gsf_crystal_eta     ;
    
    if(fabs((*gsfiter).superCluster()->eta())<1.479){//First : Barrel
      int iebhit = -1, nclust = 0;
      double amplitot = 0.0;
      double clustot = 0.0;
      
      for(reco::CaloCluster_iterator bcIt = (*gsfiter).superCluster()->clustersBegin() ; bcIt!=(*gsfiter).superCluster()->clustersEnd() ; ++bcIt) {
        // Loop over basic clusters in SC
        // bcIt seems to be a pointer to a pointer
   
        double clusterEnergy = (*bcIt)->energy();
        clustot += clusterEnergy;
        nclust += 1;
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt!=(*bcIt)->hitsAndFractions().end() ; ++rhIt) {
          // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it=EBhits->begin() ; it!=EBhits->end() ; ++it){
            // Loop over all rec hits to find the right ones
            if(rhIt->first==(*it).id()){ // Found the matching rechit
              iebhit +=1; 
              EcalRecHit hit = (*it);
              EBDetId det    = hit.id(); 
              float ampli    = hit.energy();
              amplitot += ampli;
        
              GlobalPoint poseb = geometry->getPosition(hit.detid());
              float eta_eb = poseb.eta();
              int   ieta   = det.ieta() ;      
              int   iphi   = det.iphi() ;
        
              gsf_crystal_energy  .push_back(ampli ) ;
              gsf_crystal_ietaorix.push_back(ieta  ) ;
              gsf_crystal_iphioriy.push_back(iphi  ) ;
              gsf_crystal_eta     .push_back(eta_eb) ;
            }
          }
        }
      }
    }
    //Now looking at endcaps rechits
    else{
      int ieehit = -1, nclustee = 0;
      double amplitotee = 0.0, clustotee = 0.0;
      for(reco::CaloCluster_iterator bcIt = (*gsfiter).superCluster()->clustersBegin() ; bcIt!=(*gsfiter).superCluster()->clustersEnd() ; ++bcIt){
        // Loop over basic clusters in SC
        nclustee +=1;
        double clusterEnergyee = (*bcIt)->energy();
        clustotee += clusterEnergyee;
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt!=(*bcIt)->hitsAndFractions().end() ; ++rhIt){
          // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it = EEhits->begin() ; it!=EEhits->end(); ++it){
            // Loop over all rec hits to find the right ones
            if(rhIt->first==(*it).id() ){ //found the matching rechit
              ieehit += 1;
              EcalRecHit hit = (*it);
              EEDetId det = hit.id(); 
              float ampli = hit.energy();
              amplitotee += ampli;
        
              GlobalPoint posee = geometry->getPosition(hit.detid());
              float eta_ee = posee.eta();
              int   ix     = det.ix();
              int   iy     = det.iy();

              gsf_crystal_energy  .push_back(ampli ) ;
              gsf_crystal_ietaorix.push_back(ix    ) ;
              gsf_crystal_iphioriy.push_back(iy    ) ;
              gsf_crystal_eta     .push_back(eta_ee) ;
            }
          }
        }
      }
    }
    store("HEEP_crystal_energy"  , gsf_crystal_energy  ) ;
    store("HEEP_crystal_ietaorix", gsf_crystal_ietaorix) ;
    store("HEEP_crystal_iphioriy", gsf_crystal_iphioriy) ;
    store("HEEP_crystal_eta"     , gsf_crystal_eta     ) ;
  }
  store("HEEP_nHEEP",nHeepEle);
}

void IIHEModuleHEEP::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleHEEP::beginEvent(){}
void IIHEModuleHEEP::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleHEEP::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleHEEP);
