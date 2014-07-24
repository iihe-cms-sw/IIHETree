#include "UserCode/IIHETree/interface/IIHEModuleGedGsfElectron.h"

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

IIHEModuleGedGsfElectron::IIHEModuleGedGsfElectron(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleGedGsfElectron::~IIHEModuleGedGsfElectron(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleGedGsfElectron::beginJob(){
  addBranch("gsf_n", kInt) ;
  addBranch("gsf_classification", kVectorInt) ;
  setBranchType(kVectorFloat) ;
  addBranch("gsf_energy") ;
  addBranch("gsf_p") ;
  addBranch("gsf_pt") ;
  addBranch("gsf_scE1x5") ;
  addBranch("gsf_scE5x5") ;
  addBranch("gsf_scE2x5Max") ;
  addBranch("gsf_eta") ;
  addBranch("gsf_phi") ;
  addBranch("gsf_theta") ;
  addBranch("gsf_px") ;
  addBranch("gsf_py") ;
  addBranch("gsf_pz") ;
  addBranch("gsf_deltaEtaSuperClusterTrackAtVtx") ;
  addBranch("gsf_deltaPhiSuperClusterTrackAtVtx") ;
  addBranch("gsf_hadronicOverEm") ;
  addBranch("gsf_hcalDepth1OverEcal") ;
  addBranch("gsf_hcalDepth2OverEcal") ;
  addBranch("gsf_dr03TkSumPt") ;
  addBranch("gsf_dr03EcalRecHitSumEt") ;
  addBranch("gsf_dr03HcalDepth1TowerSumEt") ;
  addBranch("gsf_dr03HcalDepth2TowerSumEt") ;
  addBranch("gsf_charge", kVectorInt) ;
  addBranch("gsf_sigmaIetaIeta") ;
  addBranch("gsf_ecaldrivenSeed"   , kVectorBool) ;
  addBranch("gsf_trackerdrivenSeed", kVectorBool) ;
  addBranch("gsf_eseffsixix") ;
  addBranch("gsf_eseffsiyiy") ;
  addBranch("gsf_eseffsirir") ;
  addBranch("gsf_preshowerEnergy") ;
  addBranch("gsf_isEB", kVectorBool) ;
  addBranch("gsf_isEE", kVectorBool) ;
  addBranch("gsf_deltaEtaSeedClusterTrackAtCalo") ;
  addBranch("gsf_deltaPhiSeedClusterTrackAtCalo") ;
  addBranch("gsf_ecalEnergy") ;
  addBranch("gsf_eSuperClusterOverP") ;
  addBranch("gsf_dxy") ;
  addBranch("gsf_dxy_beamSpot") ;
  addBranch("gsf_dxy_firstPVtx") ;
  addBranch("gsf_dxy_firstPVtxwithBS") ;
  addBranch("gsf_dxyError") ;
  addBranch("gsf_dz") ;
  addBranch("gsf_dz_beamSpot") ;
  addBranch("gsf_dz_firstPVtx") ;
  addBranch("gsf_dz_firstPVtxwithBS") ;
  addBranch("gsf_dzError") ;
  addBranch("gsf_vz") ;
  addBranch("gsf_numberOfValidHits", kVectorInt) ;
  addBranch("gsf_nLostInnerHits"   , kVectorInt) ;
  addBranch("gsf_nLostOuterHits"   , kVectorInt) ;
  addBranch("gsf_convFlags"        , kVectorInt) ;
  addBranch("gsf_convDist") ;
  addBranch("gsf_convDcot") ;
  addBranch("gsf_convRadius") ;
  addBranch("gsf_fBrem") ;
  addBranch("gsf_e1x5") ;
  addBranch("gsf_e2x5") ;
  addBranch("gsf_e5x5") ;
  addBranch("gsf_e1x3") ;
  
  setBranchType(kVectorVectorFloat) ;
  addBranch("gsf_crystal_energy"  ) ;
  addBranch("gsf_crystal_eta"     ) ;
  addBranch("gsf_eshitsixix") ;
  addBranch("gsf_eshitsiyiy") ;
  setBranchType(kVectorVectorInt) ;
  addBranch("gsf_crystal_ietaorix") ;
  addBranch("gsf_crystal_iphioriy") ;
  addBranch("gsf_hitsinfo", kVectorVectorInt) ;
  
  
  std::string branchPrefix = "gsfMatch_" ;
  std::vector<std::string> filterNames ;
  filterNames.push_back("hltL1sL1SingleEG12") ;
  filterNames.push_back("hltL1sL1Mu3p5EG12" ) ;
  filterNames.push_back("hltL1sL1SingleEG22") ;
  filterNames.push_back("hltEle33CaloIdLPixelMatchFilter") ;
  filterNames.push_back("hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter") ;
  filterNames.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter") ;
  filterNames.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") ;
  filterNames.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter") ;
  filterNames.push_back("hltEle27WP80TrackIsoFilter") ;
  filterNames.push_back("hltMu22Photon22CaloIdLHEFilter") ;
  setBranchType(kVectorBool) ;
  for(unsigned int i=0 ; i<filterNames.size() ; i++){
    std::string branchName = branchPrefix + filterNames.at(i) ;
    addBranch(branchName) ;
  }
}

// ------------ method called to for each event  ------------
void IIHEModuleGedGsfElectron::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  edm::ESHandle<CaloGeometry> pGeometry ;
  iSetup.get<CaloGeometryRecord>().get(pGeometry) ;
  CaloGeometry* geometry = (CaloGeometry*) pGeometry.product() ;
  const CaloSubdetectorGeometry* geometryES = geometry->getSubdetectorGeometry(DetId::Ecal, EcalPreshower) ;
  CaloSubdetectorTopology* topology_p = 0 ;
  if(geometryES) topology_p = new EcalPreshowerTopology(geometry) ;
  
  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
  iEvent.getByLabel("gedGsfElectrons","",pGsfElectrons);
  reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());
  
  math::XYZPoint firstpvertexwithBS(0.,0.,0.);
  Handle<reco::VertexCollection> primaryVertexCollwithBS;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS",primaryVertexCollwithBS);
  const reco::VertexCollection* pvcollwithBS = primaryVertexCollwithBS.product();

  if(pvcollwithBS->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcollwithBS->begin();
    firstpvertexwithBS.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  
  Handle<EcalRecHitCollection> EBhits;
  Handle<EcalRecHitCollection> EEhits;
  iEvent.getByLabel("reducedEcalRecHitsEB",EBhits);
  iEvent.getByLabel("reducedEcalRecHitsEE",EEhits);
  
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByLabel("offlineBeamSpot", theBeamSpot);
  math::XYZPoint beamspot(theBeamSpot->position().x(),theBeamSpot->position().y(),theBeamSpot->position().z());
  math::XYZPoint firstpvertex(0.,0.,0.);
  
  // Trigger information
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag,trigEvent);
  
  EcalClusterLazyTools lazytool(iEvent,iSetup,InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"),InputTag("reducedEcalRecHitsES"));
  
  store("gsf_n", (int)gsfelectrons.size()) ;
  for(reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin() ; gsfiter!=gsfelectrons.end() ; ++gsfiter){
    // Required for preshower variables
    reco::SuperClusterRef cl_ref = gsfiter->superCluster();
    
    //Fill the gsf related variables
    const reco::CaloClusterPtr seed = gsfiter->superCluster()->seed();
    int gsf_nLostInnerHits = gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
    int gsf_nLostOuterHits = gsfiter->gsfTrack()->trackerExpectedHitsOuter().numberOfLostHits() ;
    store("gsf_energy"                        , gsfiter->energy()) ;
    store("gsf_p"                             , gsfiter->p()) ;
    store("gsf_pt"                            , gsfiter->pt()) ;
    store("gsf_classification"                , gsfiter->classification()) ;
    store("gsf_scE1x5"                        , gsfiter->scE1x5()) ;
    store("gsf_scE5x5"                        , gsfiter->scE5x5()) ;
    store("gsf_scE2x5Max"                     , gsfiter->scE2x5Max()) ;
    store("gsf_eta"                           , gsfiter->eta()) ;
    store("gsf_phi"                           , gsfiter->phi()) ;
    store("gsf_theta"                         , gsfiter->theta());
    store("gsf_px"                            , gsfiter->px()) ;
    store("gsf_py"                            , gsfiter->py()) ;
    store("gsf_pz"                            , gsfiter->pz()) ;
    store("gsf_deltaEtaSuperClusterTrackAtVtx", gsfiter->deltaEtaSuperClusterTrackAtVtx()) ;
    store("gsf_deltaPhiSuperClusterTrackAtVtx", gsfiter->deltaPhiSuperClusterTrackAtVtx()) ;
    store("gsf_hadronicOverEm"                , gsfiter->hadronicOverEm()) ;
    store("gsf_hcalDepth1OverEcal"            , gsfiter->hcalDepth1OverEcal()) ;
    store("gsf_hcalDepth2OverEcal"            , gsfiter->hcalDepth2OverEcal()) ;
    store("gsf_dr03TkSumPt"                   , gsfiter->dr03TkSumPt()) ;
    store("gsf_dr03EcalRecHitSumEt"           , gsfiter->dr03EcalRecHitSumEt()) ;
    store("gsf_dr03HcalDepth1TowerSumEt"      , gsfiter->dr03HcalDepth1TowerSumEt()) ;
    store("gsf_dr03HcalDepth2TowerSumEt"      , gsfiter->dr03HcalDepth2TowerSumEt()) ;
    store("gsf_charge"                        , gsfiter->charge()) ;
    store("gsf_sigmaIetaIeta"                 , gsfiter->sigmaIetaIeta()) ;
    store("gsf_ecaldrivenSeed"                , gsfiter->ecalDrivenSeed()   ) ;
    store("gsf_trackerdrivenSeed"             , gsfiter->trackerDrivenSeed()) ;
    store("gsf_eseffsixix"                    , lazytool.eseffsixix(*cl_ref)) ;
    store("gsf_eseffsiyiy"                    , lazytool.eseffsiyiy(*cl_ref)) ;
    store("gsf_eseffsirir"                    , lazytool.eseffsirir(*cl_ref)) ;
    store("gsf_preshowerEnergy"               , gsfiter->superCluster()->preshowerEnergy()) ;
    store("gsf_isEB"                          , gsfiter->isEB());
    store("gsf_isEE"                          , gsfiter->isEE());
    store("gsf_deltaEtaSeedClusterTrackAtCalo", gsfiter->deltaEtaSeedClusterTrackAtCalo());
    store("gsf_deltaPhiSeedClusterTrackAtCalo", gsfiter->deltaPhiSeedClusterTrackAtCalo());
    store("gsf_ecalEnergy"                    , gsfiter->ecalEnergy());
    store("gsf_eSuperClusterOverP"            , gsfiter->eSuperClusterOverP());
    store("gsf_dxy"                           , gsfiter->gsfTrack()->dxy());
    store("gsf_dxy_beamSpot"                  , gsfiter->gsfTrack()->dxy(beamspot));
    store("gsf_dxy_firstPVtx"                 , gsfiter->gsfTrack()->dxy(firstpvertex));
    store("gsf_dxy_firstPVtxwithBS"           , gsfiter->gsfTrack()->dxy(firstpvertexwithBS));
    store("gsf_dxyError"                      , gsfiter->gsfTrack()->dxyError());
    store("gsf_dz"                            , gsfiter->gsfTrack()->dz());
    store("gsf_dz_beamSpot"                   , gsfiter->gsfTrack()->dz(beamspot));
    store("gsf_dz_firstPVtx"                  , gsfiter->gsfTrack()->dz(firstpvertex));
    store("gsf_dz_firstPVtxwithBS"            , gsfiter->gsfTrack()->dz(firstpvertexwithBS));
    store("gsf_dzError"                       , gsfiter->gsfTrack()->dzError()); 
    store("gsf_vz"                            , gsfiter->gsfTrack()->vz());
    store("gsf_numberOfValidHits"             , gsfiter->gsfTrack()->numberOfValidHits());
    store("gsf_nLostInnerHits"                , gsf_nLostInnerHits);
    store("gsf_nLostOuterHits"                , gsf_nLostOuterHits);
    store("gsf_convFlags"                     , gsfiter->convFlags());
    store("gsf_convDist"                      , gsfiter->convDist());
    store("gsf_convDcot"                      , gsfiter->convDcot());
    store("gsf_convRadius"                    , gsfiter->convRadius());
    store("gsf_fBrem"                         , gsfiter->fbrem());
    store("gsf_e1x5"                          , gsfiter->e1x5()) ;
    store("gsf_e2x5"                          , gsfiter->e2x5Max()) ;
    store("gsf_e5x5"                          , gsfiter->e5x5()) ;
    store("gsf_e1x3"                          , lazytool.e1x3(*seed)) ;
    
    //http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/TrackReco/interface/HitPattern.h?revision=1.32&view=markup
    reco::HitPattern kfHitPattern = gsfiter->gsfTrack()->hitPattern();
    int nbtrackhits = kfHitPattern.numberOfHits();
    std::vector<int> gsf_hitsinfo ;
    for(int hititer=0; hititer<25;hititer++){
      int myhitbin = (hititer<nbtrackhits) ? kfHitPattern.getHitPattern(hititer) : 0 ;
      gsf_hitsinfo.push_back(myhitbin) ;
    }
    store("gsf_hitsinfo", gsf_hitsinfo) ;
    
    // Get the preshower hits
    double x = gsfiter->superCluster()->x() ;
    double y = gsfiter->superCluster()->y() ;
    double z = gsfiter->superCluster()->z() ;
    store("gsf_eshitsixix", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 1)) ;
    store("gsf_eshitsiyiy", lazytool.getESHits(x, y, z, lazytool.rechits_map_, geometry, topology_p, 0, 2)) ;
    
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
    store("gsf_crystal_energy"  , gsf_crystal_energy  ) ;
    store("gsf_crystal_ietaorix", gsf_crystal_ietaorix) ;
    store("gsf_crystal_iphioriy", gsf_crystal_iphioriy) ;
    store("gsf_crystal_eta"     , gsf_crystal_eta     ) ;
    
    //////////////////////////////////////////////////////////////////////////////////////
    //                                 Trigger matching                                 //
    //////////////////////////////////////////////////////////////////////////////////////
    const double barrelEnd       = 1.4791;
    const double regionEtaSizeEB = 0.522 ;
    const double regionEtaSizeEE = 1.0   ;
    const double regionPhiSize   = 1.044 ;
    
    std::string branchPrefix = "gsfMatch_" ;
    std::vector<std::string> filterNames ;
    filterNames.push_back("hltL1sL1SingleEG12") ;
    filterNames.push_back("hltL1sL1Mu3p5EG12" ) ;
    filterNames.push_back("hltL1sL1SingleEG22") ;
    
    // L1 hltL1sL1SingleEG12
    // Careful that L1 triggers only have discrete eta phi. Need to be extremely loose. 
    // See here: http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/SHarper/SHNtupliser/src/SHTrigInfo.cc?revision=1.5&view=markup&pathrev=HEAD
    // It is important to specify the right HLT process for the filter, not doing this is a common bug
    for(unsigned int iFilter=0 ; iFilter<filterNames.size() ; iFilter++){
      bool gsfMatch = false ;
      int filterIndex = trigEvent->filterIndex(edm::InputTag(filterNames.at(iFilter),"",trigEventTag.process())); 
      if(filterIndex<trigEvent->sizeFilters()){ 
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
        // Now loop of the trigger objects passing filter
        for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          // Do what you want with the trigger objects, you have
          // eta,phi,pt,mass,p,px,py,pz,et,energy accessors
          
          float objeta = obj.eta(); 
          float objphi = obj.phi();
          
          double etaBinLow  = 0.0 ;
          double etaBinHigh = 0.0 ;
          
          if(fabs(objeta) < barrelEnd){
            etaBinLow  = objeta - regionEtaSizeEB/2.;
            etaBinHigh = etaBinLow + regionEtaSizeEB;
          }
          else{
            etaBinLow  = objeta - regionEtaSizeEE/2.;
            etaBinHigh = etaBinLow + regionEtaSizeEE;
          }
          
          float deltaPhi = reco::deltaPhi(gsfiter->phi(),objphi);
          
          if(gsfiter->eta() < etaBinHigh && gsfiter->eta() > etaBinLow &&   deltaPhi <regionPhiSize/2. )  {
            gsfMatch = true ;
            break ;
          }
        }
      }
      std::string branchName = branchPrefix + filterNames.at(iFilter) ;
      store(branchName, gsfMatch) ;
    }//end filter size check
    
    filterNames.clear() ;
    filterNames.push_back("hltEle33CaloIdLPixelMatchFilter") ;
    filterNames.push_back("hltDiEle33CaloIdLGsfTrkIdVLDPhiDoubleFilter") ;
    filterNames.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter") ;
    filterNames.push_back("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter") ;
    filterNames.push_back("hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter") ;
    filterNames.push_back("hltEle27WP80TrackIsoFilter") ;
    filterNames.push_back("hltMu22Photon22CaloIdLHEFilter") ;
    for(unsigned iFilter=0 ; iFilter<filterNames.size() ; iFilter++){
      bool gsfMatch = false ;
      // It is important to specify the right HLT process for the filter, not doing this is a common bug
      trigger::size_type filterIndex = trigEvent->filterIndex(edm::InputTag(filterNames.at(iFilter),"",trigEventTag.process())); 
      if(filterIndex<trigEvent->sizeFilters()){ 
        const trigger::Keys& trigKeys = trigEvent->filterKeys(filterIndex); 
        const trigger::TriggerObjectCollection & trigObjColl(trigEvent->getObjects());
        // Now loop over the trigger objects passing filter
        for(trigger::Keys::const_iterator keyIt = trigKeys.begin(); keyIt != trigKeys.end(); ++keyIt) { 
          const trigger::TriggerObject& obj = trigObjColl[*keyIt];
          if(deltaR(gsfiter->eta(),gsfiter->phi(),obj.eta(), obj.phi())<0.5){
            gsfMatch = true ;
          }
        }
      }//end filter size check
      std::string branchName = branchPrefix + filterNames.at(iFilter) ;
      store(branchName, gsfMatch) ;
    }
  }
}

void IIHEModuleGedGsfElectron::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleGedGsfElectron::beginEvent(){}
void IIHEModuleGedGsfElectron::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleGedGsfElectron::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleGedGsfElectron);
