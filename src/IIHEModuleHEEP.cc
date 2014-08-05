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
  EcalHcal1EffAreaEndcaps_ = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaEndcaps", 0.) ;
  EcalHcal1EffAreaBarrel_  = iConfig.getUntrackedParameter<double>("EcalHcal1EffAreaBarrel" , 0.) ;
}
IIHEModuleHEEP::~IIHEModuleHEEP(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleHEEP::beginJob(){
  setBranchType(kVectorBool) ;
  addBranch("HEEP_gsfpass_ET"             ) ;
  addBranch("HEEP_gsfpass_PT"             ) ;
  addBranch("HEEP_gsfpass_DETETA"         ) ;
  addBranch("HEEP_gsfpass_CRACK"          ) ;
  addBranch("HEEP_gsfpass_DETAIN"         ) ;
  addBranch("HEEP_gsfpass_DPHIIN"         ) ;
  addBranch("HEEP_gsfpass_HADEM"          ) ;
  addBranch("HEEP_gsfpass_SIGMAIETAIETA"  ) ;
  addBranch("HEEP_gsfpass_E2X5OVER5X5"    ) ;
  addBranch("HEEP_gsfpass_ISOLEMHADDEPTH1") ;
  addBranch("HEEP_gsfpass_ISOLHADDEPTH2"  ) ;
  addBranch("HEEP_gsfpass_ISOLPTTRKS"     ) ;
  addBranch("HEEP_gsfpass_ECALDRIVEN"     ) ;
  addBranch("HEEP_gsfpass_INVALID"        ) ;
  addBranch("HEEP_gsfpass_NOMISSINGHITS"  ) ;
  addBranch("HEEP_gsfpass_NOCONVERSION"   ) ;
  addBranch("HEEP_gsfpass_DXYFIRSTPV"     ) ;
  addBranch("HEEP_gsfpass_ID"             ) ;
  addBranch("HEEP_gsfpass_ISO"            ) ;
  addBranch("HEEP_gsfpass_HEEP"           ) ;
  addBranch("HEEP_nHEEP", kInt);
}

// ------------ method called to for each event  ------------
void IIHEModuleHEEP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  // rho variable
  float rho = 0 ;
  edm::Handle<double> rho_ ;
  bool isrho = iEvent.getByLabel(edm::InputTag("kt6PFJets:rho"),rho_) ;
  if(isrho) rho = *rho_ ;
  
  edm::Handle<reco::GsfElectronCollection> pGsfElectrons;
  iEvent.getByLabel("gedGsfElectrons","",pGsfElectrons);
  reco::GsfElectronCollection gsfelectrons(pGsfElectrons->begin(),pGsfElectrons->end());
  
  // Retrieve primary vertex collection
  Handle<reco::VertexCollection> primaryVertexColl;
  iEvent.getByLabel("offlinePrimaryVertices",primaryVertexColl);
  const reco::VertexCollection* pvcoll = primaryVertexColl.product();
  
  math::XYZPoint firstpvertex      (0.0,0.0,0.0) ;
  // We take only the first primary vertex, i.e. the one with the electrons
  if(pvcoll->size() > 0) {
    reco::VertexCollection::const_iterator firstpv = pvcoll->begin();
    firstpvertex.SetXYZ(firstpv->x(),firstpv->y(),firstpv->z());
  }
  
  int nHeepEle = 0 ;
  for(reco::GsfElectronCollection::const_iterator gsfiter = gsfelectrons.begin() ; gsfiter!=gsfelectrons.end() ; ++gsfiter){
    //////////////////////////////////////////////////////////////////////////////////////
    //                                    HEEP cutflow                                  //
    //////////////////////////////////////////////////////////////////////////////////////
    float gsf_gsfet           = gsfiter->caloEnergy()*sin(gsfiter->p4().theta()) ;
    float gsfsc_eta = gsfiter->superCluster()->eta() ;
    float gsf_deltaeta        = gsfiter->deltaEtaSuperClusterTrackAtVtx() ;
    float gsf_deltaphi        = gsfiter->deltaPhiSuperClusterTrackAtVtx() ;
    float gsf_hovere          = gsfiter->hadronicOverEm() ;
    float gsf_sigmaIetaIeta   = gsfiter->sigmaIetaIeta() ;
    float gsf_e2x5overe5x5    = gsfiter->scE2x5Max()/gsfiter->scE5x5() ;
    float gsf_e1x5overe5x5    = gsfiter->scE1x5()/gsfiter->scE5x5() ;
    float gsf_ecaliso         = gsfiter->dr03EcalRecHitSumEt() ;
    float gsf_hcaliso1        = gsfiter->dr03HcalDepth1TowerSumEt() ;
    float gsf_trackiso        = gsfiter->dr03TkSumPt() ;
    bool  gsf_isecaldriven    = gsfiter->ecalDrivenSeed() ;
    int   gsf_nLostInnerHits  = gsfiter->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits() ;
    int   gsf_convFlags       = gsfiter->convFlags() ;
    float gsf_dxy_firstPVtx   = gsfiter->gsfTrack()->dxy(firstpvertex) ;
    
    bool gsfetbarrel         = gsf_gsfet > 35.0 ;
    bool gsfetendcap         = gsf_gsfet > 35.0 ;
    bool barrelsc            = fabs(gsfsc_eta) < 1.442 ;
    bool endcapsc            = (fabs(gsfsc_eta) > 1.56) && (fabs(gsfsc_eta) < 2.5) ;
    bool deltaetabarrel      = fabs(gsf_deltaeta) < 0.005 ;
    bool deltaetaendcap      = fabs(gsf_deltaeta) < 0.007 ;
    bool deltaphibarrel      = fabs(gsf_deltaphi) < 0.06  ;
    bool deltaphiendcap      = fabs(gsf_deltaphi) < 0.06  ;
    bool hoverebarrel        = gsf_hovere < 0.05 ;
    bool hovereendcap        = gsf_hovere < 0.05 ;
    bool sigmaIetaIetabarrel = true ;
    bool sigmaIetaIetaendcap = gsf_sigmaIetaIeta < 0.03 ;
    bool e2x5overe5x5barrel  = (gsf_e2x5overe5x5 > 0.94) || (gsf_e1x5overe5x5 > 0.83);
    bool e2x5overe5x5endcap  = true ;
    bool ecalisobarrel       = (gsf_ecaliso+gsf_hcaliso1) < (2.+0.03*gsf_gsfet + rho*EcalHcal1EffAreaBarrel_);
    bool ecalisoendcap       = true ;
    if(gsf_gsfet<50.0){ ecalisoendcap = (gsf_ecaliso+gsf_hcaliso1) <  2.5+ rho*EcalHcal1EffAreaEndcaps_ ; }
    else              { ecalisoendcap = (gsf_ecaliso+gsf_hcaliso1) < (2.5+0.03*(gsf_gsfet-50.0)+ rho*EcalHcal1EffAreaEndcaps_ ) ; }
    bool hcaliso2barrel      = true ;
    bool hcaliso2endcap      = true ;
    bool trackisobarrel      = gsf_trackiso<5.0 ;
    bool trackisoendcap      = gsf_trackiso<5.0 ;
    bool noMissingHits       = gsf_nLostInnerHits<=1 ;
    bool noConversion        = gsf_convFlags != 3 ; 
    bool dxyfirstpvbarrel    = fabs(gsf_dxy_firstPVtx) <0.02 ;
    bool dxyfirstpvendcaps   = fabs(gsf_dxy_firstPVtx) <0.05 ;

    //Boolean HEEP cuts
    bool gsfpass_ET              = (gsfetbarrel && barrelsc) || (gsfetendcap && endcapsc) ;
    bool gsfpass_PT              = true ;
    bool gsfpass_DETETA          = true ;
    bool gsfpass_CRACK           = true ;
    bool gsfpass_DETAIN          = (deltaetabarrel && barrelsc) || (deltaetaendcap && endcapsc) ;
    bool gsfpass_DPHIIN          = (deltaphibarrel && barrelsc) || (deltaphiendcap && endcapsc) ;
    bool gsfpass_HADEM           = (hoverebarrel   && barrelsc) || (hovereendcap   && endcapsc) ;
    bool gsfpass_SIGMAIETAIETA   = (sigmaIetaIetabarrel && barrelsc) || (sigmaIetaIetaendcap && endcapsc) ;
    bool gsfpass_E2X5OVER5X5     = (e2x5overe5x5barrel  && barrelsc) || (e2x5overe5x5endcap  && endcapsc) ;
    bool gsfpass_ISOLEMHADDEPTH1 = (ecalisobarrel  && barrelsc) || (ecalisoendcap  && endcapsc) ;
    bool gsfpass_ISOLHADDEPTH2   = (hcaliso2barrel && barrelsc) || (hcaliso2endcap && endcapsc) ;
    bool gsfpass_ISOLPTTRKS      = (trackisobarrel && barrelsc) || (trackisoendcap && endcapsc) ;
    bool gsfpass_ECALDRIVEN      = gsf_isecaldriven ;
    bool gsfpass_INVALID         = true ;
    bool gsfpass_NOMISSINGHITS   = noMissingHits ;
    bool gsfpass_NOCONVERSION    = noConversion ;
    bool gsfpass_DXYFIRSTPV      = (dxyfirstpvbarrel && barrelsc) || (dxyfirstpvendcaps && endcapsc) ;
    bool gsfpass_ID              = (gsfpass_DETAIN && gsfpass_DPHIIN && gsfpass_HADEM && gsfpass_SIGMAIETAIETA && gsfpass_E2X5OVER5X5) ;
    bool gsfpass_ISO             = (gsfpass_ISOLEMHADDEPTH1 && gsfpass_ISOLHADDEPTH2 && gsfpass_ISOLPTTRKS) ;
    bool gsfpass_HEEP = gsfpass_ET && gsfpass_PT && gsfpass_DETETA && gsfpass_CRACK&& gsfpass_DETAIN && gsfpass_DPHIIN && gsfpass_HADEM && gsfpass_SIGMAIETAIETA && gsfpass_E2X5OVER5X5     && gsfpass_ISOLEMHADDEPTH1 && gsfpass_ISOLHADDEPTH2 && gsfpass_ISOLPTTRKS && gsfpass_ECALDRIVEN && gsfpass_INVALID && gsfpass_NOMISSINGHITS && gsfpass_NOCONVERSION && gsfpass_DXYFIRSTPV && gsfpass_ID && gsfpass_ISO ;
    
    store("HEEP_gsfpass_ET"             , gsfpass_ET              ) ;
    store("HEEP_gsfpass_PT"             , gsfpass_PT              ) ;
    store("HEEP_gsfpass_DETETA"         , gsfpass_DETETA          ) ;
    store("HEEP_gsfpass_CRACK"          , gsfpass_CRACK           ) ;
    store("HEEP_gsfpass_DETAIN"         , gsfpass_DETAIN          ) ;
    store("HEEP_gsfpass_DPHIIN"         , gsfpass_DPHIIN          ) ;
    store("HEEP_gsfpass_HADEM"          , gsfpass_HADEM           ) ;
    store("HEEP_gsfpass_SIGMAIETAIETA"  , gsfpass_SIGMAIETAIETA   ) ;
    store("HEEP_gsfpass_E2X5OVER5X5"    , gsfpass_E2X5OVER5X5     ) ;
    store("HEEP_gsfpass_ISOLEMHADDEPTH1", gsfpass_ISOLEMHADDEPTH1 ) ;
    store("HEEP_gsfpass_ISOLHADDEPTH2"  , gsfpass_ISOLHADDEPTH2   ) ;
    store("HEEP_gsfpass_ISOLPTTRKS"     , gsfpass_ISOLPTTRKS      ) ;
    store("HEEP_gsfpass_ECALDRIVEN"     , gsfpass_ECALDRIVEN      ) ;
    store("HEEP_gsfpass_INVALID"        , gsfpass_INVALID         ) ;
    store("HEEP_gsfpass_NOMISSINGHITS"  , gsfpass_NOMISSINGHITS   ) ;
    store("HEEP_gsfpass_NOCONVERSION"   , gsfpass_NOCONVERSION    ) ;
    store("HEEP_gsfpass_DXYFIRSTPV"     , gsfpass_DXYFIRSTPV      ) ;
    store("HEEP_gsfpass_ID"             , gsfpass_ID              ) ;
    store("HEEP_gsfpass_ISO"            , gsfpass_ISO             ) ;
    store("HEEP_gsfpass_HEEP"           , gsfpass_HEEP            ) ;
    if(gsfpass_HEEP) ++nHeepEle ;
  }
  store("HEEP_nHEEP",nHeepEle);
}

void IIHEModuleHEEP::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleHEEP::beginEvent(){}
void IIHEModuleHEEP::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleHEEP::endJob(){}

DEFINE_FWK_MODULE(IIHEModuleHEEP);
