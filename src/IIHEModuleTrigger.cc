#include "UserCode/IIHETree/interface/IIHEModuleTrigger.h"

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

IIHEModuleTrigger::IIHEModuleTrigger(const edm::ParameterSet& iConfig): IIHEModule(iConfig){}
IIHEModuleTrigger::~IIHEModuleTrigger(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleTrigger::beginJob(){
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
void IIHEModuleTrigger::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  
  // Trigger information
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","HLT");
  edm::Handle<trigger::TriggerEvent> trigEvent ; 
  iEvent.getByLabel(trigEventTag,trigEvent) ;
  
  for(reco::GsfElectronCollection::const_iterator gsfiter=electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
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
            etaBinLow  = objeta    - regionEtaSizeEB/2.0 ;
            etaBinHigh = etaBinLow + regionEtaSizeEB     ;
          }
          else{
            etaBinLow  = objeta    - regionEtaSizeEE/2.0 ;
            etaBinHigh = etaBinLow + regionEtaSizeEE     ;
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
    }
    
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

void IIHEModuleTrigger::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleTrigger::beginEvent(){}
void IIHEModuleTrigger::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleTrigger::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleTrigger);
