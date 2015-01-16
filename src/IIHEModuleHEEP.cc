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
  EcalHcal1EffAreaBarrel_  = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaBarrel" , 0.28) ;
  EcalHcal1EffAreaEndcaps_ = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaEndcaps", 0.28) ;
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
  //                                          4.1                                       //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_41 = "HEEP_cutflow41" ;
  HEEPCutflow_41_ID_        = new HEEPCutCollection(prefix_41 + "_ID"       , this) ;
  HEEPCutflow_41_isolation_ = new HEEPCutCollection(prefix_41 + "_isolation", this) ;
  HEEPCutflow_41_total_     = new HEEPCutCollection(prefix_41 + "_total"    , this) ;
  
  // These cuts must be updated so declare them separately
  cut_41_isolEMHadDepth1_       = new HEEPCut_41_isolEMHadDepth1(prefix_41 + "_isolEMHadDepth1", this) ;
  cut_41_dxyFirstPV_            = new HEEPCut_41_dxyFirstPV      (prefix_41 + "_dxyFirstPV"    , this) ;
  
  cut_41_isolEMHadDepth1_->setEcalHcal1EffAreaBarrel (EcalHcal1EffAreaBarrel_ ) ;
  cut_41_isolEMHadDepth1_->setEcalHcal1EffAreaEndcaps(EcalHcal1EffAreaEndcaps_) ;
  cut_41_isolEMHadDepth1_->setRho(rho_) ;
  
  // Define the ID
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_Et           (prefix_41 + "_Et"           , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_eta          (prefix_41 + "_eta"          , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_EcalDriven   (prefix_41 + "_EcalDriven"   , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_dEtaIn       (prefix_41 + "_dEtaIn"       , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_dPhiIn       (prefix_41 + "_dPhiIn"       , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_HOverE       (prefix_41 + "_HOverE"       , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_SigmaIetaIeta(prefix_41 + "_SigmaIetaIeta", this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_E1x5OverE5x5 (prefix_41 + "_E1x5OverE5x5" , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_E2x5OverE5x5 (prefix_41 + "_E2x5OverE5x5" , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_missingHits  (prefix_41 + "_missingHits"  , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) cut_41_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_41_isolation_->addCut( (HEEPCutBase*) cut_41_isolEMHadDepth1_) ;
  HEEPCutflow_41_isolation_->addCut( (HEEPCutBase*) new HEEPCut_41_IsolPtTrks(prefix_41 + "_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_41_total_->addCutCollection(HEEPCutflow_41_ID_       ) ;
  HEEPCutflow_41_total_->addCutCollection(HEEPCutflow_41_isolation_) ;
  
  HEEPCutflow_41_total_->config() ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                        5.0 50ns                                    //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_50_50ns = "HEEP_cutflow50_50ns" ;
  HEEPCutflow_50_50ns_ID_        = new HEEPCutCollection(prefix_50_50ns+"_ID"       , this) ;
  HEEPCutflow_50_50ns_isolation_ = new HEEPCutCollection(prefix_50_50ns+"_isolation", this) ;
  HEEPCutflow_50_50ns_total_     = new HEEPCutCollection(prefix_50_50ns+"_total"    , this) ;
  
  // These cuts must be updated so declare them separately
  cut_50_50ns_isolEMHadDepth1_ = new HEEPCut_50_50ns_isolEMHadDepth1(prefix_50_50ns+"_isolEMHadDepth1", this) ;
  cut_50_50ns_dxyFirstPV_      = new HEEPCut_50_50ns_dxyFirstPV     (prefix_50_50ns+"_dxyFirstPV"     , this) ;
  
  cut_50_50ns_isolEMHadDepth1_->setEcalHcal1EffAreaBarrel (EcalHcal1EffAreaBarrel_ ) ;
  cut_50_50ns_isolEMHadDepth1_->setEcalHcal1EffAreaEndcaps(EcalHcal1EffAreaEndcaps_) ;
  cut_50_50ns_isolEMHadDepth1_->setRho(rho_) ;
  
  // Define the ID
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_Et           (prefix_50_50ns+"_Et"           , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_eta          (prefix_50_50ns+"_eta"          , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_EcalDriven   (prefix_50_50ns+"_EcalDriven"   , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_dEtaIn       (prefix_50_50ns+"_dEtaIn"       , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_dPhiIn       (prefix_50_50ns+"_dPhiIn"       , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_HOverE       (prefix_50_50ns+"_HOverE"       , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_SigmaIetaIeta(prefix_50_50ns+"_SigmaIetaIeta", this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_E1x5OverE5x5 (prefix_50_50ns+"_E1x5OverE5x5" , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_E2x5OverE5x5 (prefix_50_50ns+"_E2x5OverE5x5" , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_missingHits  (prefix_50_50ns+"_missingHits"  , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) cut_50_50ns_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_50_50ns_isolation_->addCut( (HEEPCutBase*) cut_50_50ns_isolEMHadDepth1_) ;
  HEEPCutflow_50_50ns_isolation_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_IsolPtTrks(prefix_50_50ns+"_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_50_50ns_total_->addCutCollection(HEEPCutflow_50_50ns_ID_       ) ;
  HEEPCutflow_50_50ns_total_->addCutCollection(HEEPCutflow_50_50ns_isolation_) ;
  
  HEEPCutflow_50_50ns_total_->config() ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                        5.0 25ns                                    //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_50_25ns = "HEEP_cutflow50_25ns" ;
  HEEPCutflow_50_25ns_ID_        = new HEEPCutCollection(prefix_50_25ns+"_ID"       , this) ;
  HEEPCutflow_50_25ns_isolation_ = new HEEPCutCollection(prefix_50_25ns+"_isolation", this) ;
  HEEPCutflow_50_25ns_total_     = new HEEPCutCollection(prefix_50_25ns+"_total"    , this) ;
  
  // These cuts must be updated so declare them separately
  cut_50_25ns_isolEMHadDepth1_ = new HEEPCut_50_25ns_isolEMHadDepth1(prefix_50_25ns+"_isolEMHadDepth1", this) ;
  cut_50_25ns_dxyFirstPV_      = new HEEPCut_50_25ns_dxyFirstPV     (prefix_50_25ns+"_dxyFirstPV"     , this) ;
  
  cut_50_25ns_isolEMHadDepth1_->setEcalHcal1EffAreaBarrel (EcalHcal1EffAreaBarrel_ ) ;
  cut_50_25ns_isolEMHadDepth1_->setEcalHcal1EffAreaEndcaps(EcalHcal1EffAreaEndcaps_) ;
  cut_50_25ns_isolEMHadDepth1_->setRho(rho_) ;
  
  // Define the ID
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_Et           (prefix_50_25ns+"_Et"           , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_eta          (prefix_50_25ns+"_eta"          , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_EcalDriven   (prefix_50_25ns+"_EcalDriven"   , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_dEtaIn       (prefix_50_25ns+"_dEtaIn"       , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_dPhiIn       (prefix_50_25ns+"_dPhiIn"       , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_HOverE       (prefix_50_25ns+"_HOverE"       , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_SigmaIetaIeta(prefix_50_25ns+"_SigmaIetaIeta", this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_E1x5OverE5x5 (prefix_50_25ns+"_E1x5OverE5x5" , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_E2x5OverE5x5 (prefix_50_25ns+"_E2x5OverE5x5" , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_missingHits  (prefix_50_25ns+"_missingHits"  , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) cut_50_25ns_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_50_25ns_isolation_->addCut( (HEEPCutBase*) cut_50_25ns_isolEMHadDepth1_) ;
  HEEPCutflow_50_25ns_isolation_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_IsolPtTrks(prefix_50_25ns+"_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_50_25ns_total_->addCutCollection(HEEPCutflow_50_25ns_ID_       ) ;
  HEEPCutflow_50_25ns_total_->addCutCollection(HEEPCutflow_50_25ns_isolation_) ;
  
  HEEPCutflow_50_25ns_total_->config() ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                     Add everything to the vector of cutflows                       //
  ////////////////////////////////////////////////////////////////////////////////////////
  
  HEEPCutflows_.push_back(HEEPCutflow_41_ID_            ) ;
  HEEPCutflows_.push_back(HEEPCutflow_41_isolation_     ) ;
  HEEPCutflows_.push_back(HEEPCutflow_41_total_         ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_25ns_ID_       ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_25ns_isolation_) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_25ns_total_    ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_50ns_ID_       ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_50ns_isolation_) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_50ns_total_    ) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                      Triggers                                      //
  ////////////////////////////////////////////////////////////////////////////////////////
  addTriggerL1Electron("hltL1sL1SingleEG12") ;
  addTriggerL1Electron("hltL1sL1Mu3p5EG12" ) ;
  addTriggerL1Electron("hltL1sL1SingleEG22") ;
  
  float DeltaRCut = 0.5 ;
  addTriggerHLTElectron("hltEle33CaloIdLPixelMatchFilter"                               , DeltaRCut) ;
  addTriggerHLTElectron("hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter"                   , DeltaRCut) ;
  addTriggerHLTElectron("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter"      , DeltaRCut) ;
  addTriggerHLTElectron("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter", DeltaRCut) ;
  addTriggerHLTElectron("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter"        , DeltaRCut) ;
  addTriggerHLTElectron("hltEle27WP80TrackIsoFilter"                                    , DeltaRCut) ;
  addTriggerHLTElectron("hltMu22Photon22CaloIdLHEFilter"                                , DeltaRCut) ;
}

// ------------ method called to for each event  ------------
void IIHEModuleHEEP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  
  // Pass parameters to the HEEP cutflow objects
  math::XYZPoint* firstPrimaryVertex = parent_->getFirstPrimaryVertex() ;
  cut_41_dxyFirstPV_     ->setFirstPV(firstPrimaryVertex) ;
  cut_50_50ns_dxyFirstPV_->setFirstPV(firstPrimaryVertex) ;
  cut_50_25ns_dxyFirstPV_->setFirstPV(firstPrimaryVertex) ;
  
  // Get the hit information
  Handle<EcalRecHitCollection> EBhits;
  Handle<EcalRecHitCollection> EEhits;
  iEvent.getByLabel("reducedEcalRecHitsEB",EBhits);
  iEvent.getByLabel("reducedEcalRecHitsEE",EEhits); 
  
  // Information for preshower
  edm::ESHandle<CaloGeometry> pGeometry ;
  iSetup.get<CaloGeometryRecord>().get(pGeometry) ;
  CaloGeometry* geometry = (CaloGeometry*) pGeometry.product() ;
  const CaloSubdetectorGeometry* geometryES = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower) ;
  CaloSubdetectorTopology* topology_p = 0 ;
  if(geometryES) topology_p = new EcalPreshowerTopology(geometry) ;
  
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  
  EcalClusterLazyTools lazytool(iEvent,iSetup, parent_->getReducedBarrelRecHitCollectionToken(), parent_->getReducedEndcapRecHitCollectionToken()) ;
  
  for(reco::GsfElectronCollection::const_iterator gsfiter = electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
    reco::GsfElectron* gsf = (reco::GsfElectron*) &*gsfiter ;
    
    // Required for preshower variables
    reco::SuperClusterRef    cl_ref = gsf->superCluster() ;
    const reco::CaloClusterPtr seed = gsf->superCluster()->seed() ;
    
    HEEPCutflow_41_total_     ->applyCuts(gsf) ;
    HEEPCutflow_50_50ns_total_->applyCuts(gsf) ;
    HEEPCutflow_50_25ns_total_->applyCuts(gsf) ;
    
    // Preshower information
    // Get the preshower hits
    double x = gsf->superCluster()->x() ;
    double y = gsf->superCluster()->y() ;
    double z = gsf->superCluster()->z() ;
    store("HEEP_eshitsixix", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 1)) ;
    store("HEEP_eshitsiyiy", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 2)) ;
    store("HEEP_preshowerEnergy", gsf->superCluster()->preshowerEnergy()) ;
    
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
    
    if(fabs((*gsfiter).superCluster()->eta())<1.479){ //First : Barrel
      for(reco::CaloCluster_iterator bcIt = gsf->superCluster()->clustersBegin() ; bcIt!=gsf->superCluster()->clustersEnd() ; ++bcIt) {
        // Loop over basic clusters in SC
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt!=(*bcIt)->hitsAndFractions().end() ; ++rhIt){
          // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it=EBhits->begin() ; it!=EBhits->end() ; ++it){
            // Loop over all rec hits to find the right ones
            if(rhIt->first==(*it).id()){ // Found the matching rechit
              EBDetId det    = it->id(); 
              float ampli    = it->energy();
        
              GlobalPoint poseb = geometry->getPosition(it->detid());
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
    else{ // Now looking at endcaps rechits
      for(reco::CaloCluster_iterator bcIt = gsf->superCluster()->clustersBegin() ; bcIt!=gsf->superCluster()->clustersEnd() ; ++bcIt){
        // Loop over basic clusters in SC
        for(std::vector< std::pair<DetId, float> >::const_iterator rhIt = (*bcIt)->hitsAndFractions().begin() ; rhIt!=(*bcIt)->hitsAndFractions().end() ; ++rhIt){
          // Loop over rec hits in basic cluster
          for(EcalRecHitCollection::const_iterator it = EEhits->begin() ; it!=EEhits->end(); ++it){
            // Loop over all rec hits to find the right ones
            if(rhIt->first==it->id()){ //found the matching rechit
              EEDetId det = it->id(); 
              float ampli = it->energy();
        
              GlobalPoint posee = geometry->getPosition(it->detid());
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
}

void IIHEModuleHEEP::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleHEEP::beginEvent(){
  for(unsigned i=0 ; i<HEEPCutflows_.size() ; i++){
    HEEPCutflows_.at(i)->beginEvent() ;
  }
}
void IIHEModuleHEEP::endEvent(){
  for(unsigned i=0 ; i<HEEPCutflows_.size() ; i++){
    HEEPCutflows_.at(i)->endEvent() ;
  }
}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleHEEP::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleHEEP);
