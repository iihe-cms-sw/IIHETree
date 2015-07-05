#include "UserCode/IIHETree/interface/IIHEModuleHEEP.h"
#include "UserCode/IIHETree/interface/HEEPCut.h"

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

void IIHEModuleHEEP::addHEEPParameter(int type, std::string name, std::string pset_name, float defaultValue){
  parameters_.push_back(new HEEPParameter(type, name, pset_name, defaultValue)) ;
}

IIHEModuleHEEP::IIHEModuleHEEP(const edm::ParameterSet& iConfig): IIHEModule(iConfig){
  // Decide whether to store the cutflow variables for each cutflow
  storeHEEP41_    = iConfig.getUntrackedParameter<bool>("storeHEEP41"    , true ) ;
  storeHEEP50_50_ = iConfig.getUntrackedParameter<bool>("storeHEEP50_50" , true ) ;
  storeHEEP50_25_ = iConfig.getUntrackedParameter<bool>("storeHEEP50_25" , true ) ;
  storeHEEP50_    = iConfig.getUntrackedParameter<bool>("storeHEEP50"    , true ) ;
  storeHEEP51_    = iConfig.getUntrackedParameter<bool>("storeHEEP51"    , true ) ;
  storeHEEP60_    = iConfig.getUntrackedParameter<bool>("storeHEEP60"    , true ) ;
  
  nAccept_ = 0 ;
  
  // Acceptance
  addHEEPParameter(kHC4      , "EtThresholdBarrel", "HEEP_EtThresholdBarrel_4"      , 35.0  ) ;
  addHEEPParameter(kHC41     , "EtThresholdBarrel", "HEEP_EtThresholdBarrel_41"     , 35.0  ) ;
  addHEEPParameter(kHC5      , "EtThresholdBarrel", "HEEP_EtThresholdBarrel_5"      , 35.0  ) ;
  addHEEPParameter(kHC50_50ns, "EtThresholdBarrel", "HEEP_EtThresholdBarrel_50_50ns", 35.0  ) ;
  addHEEPParameter(kHC50_25ns, "EtThresholdBarrel", "HEEP_EtThresholdBarrel_50_25ns", 35.0  ) ;
  addHEEPParameter(kHC50     , "EtThresholdBarrel", "HEEP_EtThresholdBarrel_50"     , 35.0  ) ;
  addHEEPParameter(kHC51     , "EtThresholdBarrel", "HEEP_EtThresholdBarrel_51"     , 35.0  ) ;
  addHEEPParameter(kHC60     , "EtThresholdBarrel", "HEEP_EtThresholdBarrel_60"     , 35.0  ) ;
  
  addHEEPParameter(kHC4      , "EtThresholdEndcap", "HEEP_EtThresholdBarrel_4"      , 35.0  ) ;
  addHEEPParameter(kHC41     , "EtThresholdEndcap", "HEEP_EtThresholdBarrel_41"     , 35.0  ) ;
  addHEEPParameter(kHC5      , "EtThresholdEndcap", "HEEP_EtThresholdBarrel_5"      , 35.0  ) ;
  addHEEPParameter(kHC50_50ns, "EtThresholdEndcap", "HEEP_EtThresholdBarrel_50_50ns", 35.0  ) ;
  addHEEPParameter(kHC50_25ns, "EtThresholdEndcap", "HEEP_EtThresholdBarrel_50_25ns", 35.0  ) ;
  addHEEPParameter(kHC50     , "EtThresholdEndcap", "HEEP_EtThresholdBarrel_50"     , 35.0  ) ;
  addHEEPParameter(kHC51     , "EtThresholdEndcap", "HEEP_EtThresholdBarrel_51"     , 35.0  ) ;
  addHEEPParameter(kHC60     , "EtThresholdEndcap", "HEEP_EtThresholdBarrel_60"     , 35.0  ) ;
  
  addHEEPParameter(kHC4      , "barrelEtaUpper"   , "HEEP_barrelEtaUpper_4"         , 1.442 ) ;
  addHEEPParameter(kHC41     , "barrelEtaUpper"   , "HEEP_barrelEtaUpper_41"        , 1.442 ) ;
  addHEEPParameter(kHC5      , "barrelEtaUpper"   , "HEEP_barrelEtaUpper_5"         , 1.4442) ;
  addHEEPParameter(kHC50_50ns, "barrelEtaUpper"   , "HEEP_barrelEtaUpper_50_50ns"   , 1.4442) ;
  addHEEPParameter(kHC50_25ns, "barrelEtaUpper"   , "HEEP_barrelEtaUpper_50_25ns"   , 1.4442) ;
  addHEEPParameter(kHC50     , "barrelEtaUpper"   , "HEEP_barrelEtaUpper_50"        , 1.4442) ;
  addHEEPParameter(kHC51     , "barrelEtaUpper"   , "HEEP_barrelEtaUpper_51"        , 1.4442) ;
  addHEEPParameter(kHC60     , "barrelEtaUpper"   , "HEEP_barrelEtaUpper_60"        , 1.4442) ;
  
  addHEEPParameter(kHC4      , "endcapEtaLower"   , "HEEP_endcapEtaLower_4"         , 1.56  ) ;
  addHEEPParameter(kHC41     , "endcapEtaLower"   , "HEEP_endcapEtaLower_41"        , 1.56  ) ;
  addHEEPParameter(kHC5      , "endcapEtaLower"   , "HEEP_endcapEtaLower_5"         , 1.556 ) ;
  addHEEPParameter(kHC50_50ns, "endcapEtaLower"   , "HEEP_endcapEtaLower_50_50ns"   , 1.556 ) ;
  addHEEPParameter(kHC50_25ns, "endcapEtaLower"   , "HEEP_endcapEtaLower_50_25ns"   , 1.556 ) ;
  addHEEPParameter(kHC50     , "endcapEtaLower"   , "HEEP_endcapEtaLower_50"        , 1.556 ) ;
  addHEEPParameter(kHC51     , "endcapEtaLower"   , "HEEP_endcapEtaLower_51"        , 1.556 ) ;
  addHEEPParameter(kHC60     , "endcapEtaLower"   , "HEEP_endcapEtaLower_60"        , 1.556 ) ;
  
  addHEEPParameter(kHC4      , "endcapEtaUpper"   , "HEEP_endcapEtaUpper_4"         , 2.5   ) ;
  addHEEPParameter(kHC41     , "endcapEtaUpper"   , "HEEP_endcapEtaUpper_41"        , 2.5   ) ;
  addHEEPParameter(kHC5      , "endcapEtaUpper"   , "HEEP_endcapEtaUpper_5"         , 2.5   ) ;
  addHEEPParameter(kHC50_50ns, "endcapEtaUpper"   , "HEEP_endcapEtaUpper_50_250ns"  , 2.5   ) ;
  addHEEPParameter(kHC50_25ns, "endcapEtaUpper"   , "HEEP_endcapEtaUpper_50_25ns"   , 2.5   ) ;
  addHEEPParameter(kHC50     , "endcapEtaUpper"   , "HEEP_endcapEtaUpper_50"        , 2.5   ) ;
  addHEEPParameter(kHC51     , "endcapEtaUpper"   , "HEEP_endcapEtaUpper_51"        , 2.5   ) ;
  addHEEPParameter(kHC60     , "endcapEtaUpper"   , "HEEP_endcapEtaUpper_60"        , 2.5   ) ;
  
  // ID
  addHEEPParameter(kHC41     , "dEtaInThresholdBarrel"                  , "HEEP_dEtaInThresholdBarrel_41"                        , 0.005  ) ;
  addHEEPParameter(kHC41     , "dEtaInThresholdEndcap"                  , "HEEP_dEtaInThresholdEndcap_41"                        , 0.007  ) ;
  addHEEPParameter(kHC50_50ns, "dEtaInConstantTermBarrel"               , "HEEP_dEtaInConstantTermBarrel_50_50ns"                , 0.016  ) ;
  addHEEPParameter(kHC50_50ns, "dEtaInLinearTermBarrel"                 , "HEEP_dEtaInLinearTermBarrel_50_50ns"                  , 0.0001 ) ;
  addHEEPParameter(kHC50_50ns, "dEtaInCutoffTermBarrel"                 , "HEEP_dEtaInCutoffTermBarrel_50_50ns"                  , 0.004  ) ;
  addHEEPParameter(kHC50_50ns, "dEtaInThresholdEndcap"                  , "HEEP_dEtaInThresholdEndcap_50_50ns"                   , 0.02   ) ;
  addHEEPParameter(kHC50_25ns, "dEtaInConstantTermBarrel"               , "HEEP_dEtaInConstantTermBarrel_50_25ns"                , 0.016  ) ;
  addHEEPParameter(kHC50_25ns, "dEtaInLinearTermBarrel"                 , "HEEP_dEtaInLinearTermBarrel_50_25ns"                  , 0.0001 ) ;
  addHEEPParameter(kHC50_25ns, "dEtaInCutoffTermBarrel"                 , "HEEP_dEtaInCutoffTermBarrel_50_25ns"                  , 0.004  ) ;
  addHEEPParameter(kHC50_25ns, "dEtaInConstantTermEndcap"               , "HEEP_dEtaInConstantTermEndcap_50_25ns"                , 0.016  ) ;
  addHEEPParameter(kHC50_25ns, "dEtaInLinearTermEndcap"                 , "HEEP_dEtaInLinearTermEndcap_50_25ns"                  , 0.0001 ) ;
  addHEEPParameter(kHC50_25ns, "dEtaInCutoffTermEndcap"                 , "HEEP_dEtaInCutoffTermEndcap_50_25ns"                  , 0.004  ) ;
  addHEEPParameter(kHC50     , "dEtaInConstantTermBarrel"               , "HEEP_dEtaInConstantTermBarrel_50"                     , 0.016  ) ;
  addHEEPParameter(kHC50     , "dEtaInLinearTermBarrel"                 , "HEEP_dEtaInLinearTermBarrel_50"                       , 0.0001 ) ;
  addHEEPParameter(kHC50     , "dEtaInCutoffTermBarrel"                 , "HEEP_dEtaInCutoffTermBarrel_50"                       , 0.004  ) ;
  addHEEPParameter(kHC50     , "dEtaInConstantTermEndcap"               , "HEEP_dEtaInConstantTermEndcap_50"                     , 0.015  ) ;
  addHEEPParameter(kHC50     , "dEtaInLinearTermEndcap"                 , "HEEP_dEtaInLinearTermEndcap_50"                       , 0.00085) ;
  addHEEPParameter(kHC50     , "dEtaInCutoffTermEndcap"                 , "HEEP_dEtaInCutoffTermEndcap_50"                       , 0.006  ) ;
  addHEEPParameter(kHC51     , "dEtaInThresholdBarrel"                  , "HEEP_dEtaInThresholdBarrel_51"                        , 0.004  ) ;
  addHEEPParameter(kHC51     , "dEtaInThresholdEndcap"                  , "HEEP_dEtaInThresholdEndcap_51"                        , 0.006  ) ;
  addHEEPParameter(kHC60     , "dEtaInThresholdBarrel"                  , "HEEP_dEtaInThresholdBarrel_60"                        , 0.004  ) ;
  addHEEPParameter(kHC60     , "dEtaInThresholdEndcap"                  , "HEEP_dEtaInThresholdEndcap_60"                        , 0.006  ) ;
  
  addHEEPParameter(kHC41     , "dPhiInThresholdBarrel"                  , "HEEP_dPhiInThresholdBarrel_41"                        , 0.06) ;
  addHEEPParameter(kHC41     , "dPhiInThresholdEndcap"                  , "HEEP_dPhiInThresholdEndcap_41"                        , 0.06) ;
  addHEEPParameter(kHC50_50ns, "dPhiInThresholdBarrel"                  , "HEEP_dPhiInThresholdBarrel_50_50ns"                   , 0.06) ;
  addHEEPParameter(kHC50_50ns, "dPhiInThresholdEndcap"                  , "HEEP_dPhiInThresholdEndcap_50_50ns"                   , 0.15) ;
  addHEEPParameter(kHC50_25ns, "dPhiInThresholdBarrel"                  , "HEEP_dPhiInThresholdBarrel_50_25ns"                   , 0.06) ;
  addHEEPParameter(kHC50_25ns, "dPhiInThresholdEndcap"                  , "HEEP_dPhiInThresholdEndcap_50_25ns"                   , 0.06) ;
  addHEEPParameter(kHC50     , "dPhiInThresholdBarrel"                  , "HEEP_dPhiInThresholdBarrel_50"                        , 0.06) ;
  addHEEPParameter(kHC50     , "dPhiInThresholdEndcap"                  , "HEEP_dPhiInThresholdEndcap_50"                        , 0.06) ;
  addHEEPParameter(kHC51     , "dPhiInThresholdBarrel"                  , "HEEP_dPhiInThresholdBarrel_51"                        , 0.06) ;
  addHEEPParameter(kHC51     , "dPhiInThresholdEndcap"                  , "HEEP_dPhiInThresholdEndcap_51"                        , 0.06) ;
  addHEEPParameter(kHC60     , "dPhiInThresholdBarrel"                  , "HEEP_dPhiInThresholdBarrel_60"                        , 0.06) ;
  addHEEPParameter(kHC60     , "dPhiInThresholdEndcap"                  , "HEEP_dPhiInThresholdEndcap_60"                        , 0.06) ;
  
  addHEEPParameter(kHC41     , "HOverEThresholdBarrel"                  , "HEEP_HOverEThresholdBarrel_41"                        , 0.05) ;
  addHEEPParameter(kHC41     , "HOverEThresholdEndcap"                  , "HEEP_HOverEThresholdEndcap_41"                        , 0.05) ;
  addHEEPParameter(kHC50_50ns, "HOverEReciprocalTermBarrel"             , "HEEP_HOverEReciprocalTermBarrel_50_50ns"              , 2.0 ) ;
  addHEEPParameter(kHC50_50ns, "HOverEConstantTermBarrel"               , "HEEP_HOverEConstantTermBarrel_50_50ns"                , 0.05) ;
  addHEEPParameter(kHC50_50ns, "HOverEReciprocalTermEndcap"             , "HEEP_HOverEReciprocalTermEndcap_50_50ns"              , 12.5) ;
  addHEEPParameter(kHC50_50ns, "HOverEConstantTermEndcap"               , "HEEP_HOverEConstantTermEndcap_50_50ns"                , 0.05) ;
  addHEEPParameter(kHC50_25ns, "HOverEReciprocalTermBarrel"             , "HEEP_HOverEReciprocalTermBarrel_50_25ns"              , 2.0 ) ;
  addHEEPParameter(kHC50_25ns, "HOverEConstantTermBarrel"               , "HEEP_HOverEConstantTermBarrel_50_25ns"                , 0.05) ;
  addHEEPParameter(kHC50_25ns, "HOverEReciprocalTermEndcap"             , "HEEP_HOverEReciprocalTermEndcap_50_25ns"              , 12.5) ;
  addHEEPParameter(kHC50_25ns, "HOverEConstantTermEndcap_"              , "HEEP_HOverEConstantTermEndcap_50_25ns"                , 0.05) ;
  addHEEPParameter(kHC50     , "HOverEReciprocalTermBarrel"             , "HEEP_HOverEReciprocalTermBarrel_50"                   , 2.0 ) ;
  addHEEPParameter(kHC50     , "HOverEConstantTermBarrel"               , "HEEP_HOverEConstantTermBarrel_50"                     , 0.05) ;
  addHEEPParameter(kHC50     , "HOverEReciprocalTermEndcap"             , "HEEP_HOverEReciprocalTermEndcap_50"                   , 12.5) ;
  addHEEPParameter(kHC50     , "HOverEConstantTermEndcap"               , "HEEP_HOverEConstantTermEndcap_50"                     , 0.05) ;
  addHEEPParameter(kHC51     , "HOverEReciprocalTermBarrel"             , "HEEP_HOverEReciprocalTermBarrel_51"                   , 2.0 ) ;
  addHEEPParameter(kHC51     , "HOverEConstantTermBarrel"               , "HEEP_HOverEConstantTermBarrel_51"                     , 0.05) ;
  addHEEPParameter(kHC51     , "HOverEReciprocalTermEndcap"             , "HEEP_HOverEReciprocalTermEndcap_51"                   , 12.5) ;
  addHEEPParameter(kHC51     , "HOverEConstantTermEndcap_"              , "HEEP_HOverEConstantTermEndcap_51"                     , 0.05) ;
  addHEEPParameter(kHC60     , "HOverEReciprocalTermBarrel"             , "HEEP_HOverEReciprocalTermBarrel_60"                   , 1.0 ) ;
  addHEEPParameter(kHC60     , "HOverEConstantTermBarrel"               , "HEEP_HOverEConstantTermBarrel_60"                     , 0.05) ;
  addHEEPParameter(kHC60     , "HOverEReciprocalTermEndcap"             , "HEEP_HOverEReciprocalTermEndcap_60"                   , 5.0 ) ;
  addHEEPParameter(kHC60     , "HOverEConstantTermEndcap_"              , "HEEP_HOverEConstantTermEndcap_60"                     , 0.05) ;
  
  addHEEPParameter(kHC41     , "SigmaIetaIetaThreshold"                 , "HEEP_SigmaIetaIetaThreshold_41"                       , 0.03) ;
  addHEEPParameter(kHC50_50ns, "SigmaIetaIetaThreshold"                 , "HEEP_SigmaIetaIetaThreshold_50_50ns"                  , 0.03) ;
  addHEEPParameter(kHC50_25ns, "SigmaIetaIetaThreshold"                 , "HEEP_SigmaIetaIetaThreshold_50_25ns"                  , 0.03) ;
  addHEEPParameter(kHC50     , "SigmaIetaIetaThreshold"                 , "HEEP_SigmaIetaIetaThreshold_50"                       , 0.03) ;
  addHEEPParameter(kHC51     , "SigmaIetaIetaThreshold"                 , "HEEP_SigmaIetaIetaThreshold_50"                       , 0.03) ;
  addHEEPParameter(kHC60     , "SigmaIetaIetaThreshold"                 , "HEEP_SigmaIetaIetaThreshold_60"                       , 0.03) ;
  
  addHEEPParameter(kHC41     , "E1x5threshold"                          , "HEEP_E1x5threshold_41"                                , 0.83) ;
  addHEEPParameter(kHC41     , "E2x5threshold"                          , "HEEP_E2x5threshold_41"                                , 0.94) ;
  addHEEPParameter(kHC50_50ns, "E1x5threshold"                          , "HEEP_E1x5threshold_50_50ns"                           , 0.83) ;
  addHEEPParameter(kHC50_50ns, "E2x5threshold"                          , "HEEP_E2x5threshold_50_50ns"                           , 0.94) ;
  addHEEPParameter(kHC50_25ns, "E1x5threshold"                          , "HEEP_E1x5threshold_50_25ns"                           , 0.83) ;
  addHEEPParameter(kHC50_25ns, "E2x5threshold"                          , "HEEP_E2x5threshold_50_25ns"                           , 0.94) ;
  addHEEPParameter(kHC50     , "E1x5threshold"                          , "HEEP_E1x5threshold_50"                                , 0.83) ;
  addHEEPParameter(kHC50     , "E2x5threshold"                          , "HEEP_E2x5threshold_50"                                , 0.94) ;
  addHEEPParameter(kHC51     , "E1x5threshold"                          , "HEEP_E1x5threshold_51"                                , 0.83) ;
  addHEEPParameter(kHC51     , "E2x5threshold"                          , "HEEP_E2x5threshold_51"                                , 0.94) ;
  addHEEPParameter(kHC60     , "E1x5threshold"                          , "HEEP_E1x5threshold_60"                                , 0.83) ;
  addHEEPParameter(kHC60     , "E2x5threshold"                          , "HEEP_E2x5threshold_60"                                , 0.94) ;
  
  addHEEPParameter(kHC41     , "missingHitsThreshold"                   , "HEEP_missingHitsThreshold_41"                         , 1   ) ;
  addHEEPParameter(kHC50_50ns, "missingHitsThreshold"                   , "HEEP_missingHitsThreshold_50_50ns"                    , 1   ) ;
  addHEEPParameter(kHC50_25ns, "missingHitsThreshold"                   , "HEEP_missingHitsThreshold_50_25ns"                    , 1   ) ;
  addHEEPParameter(kHC50     , "missingHitsThreshold"                   , "HEEP_missingHitsThreshold_50"                         , 1   ) ;
  addHEEPParameter(kHC51     , "missingHitsThreshold"                   , "HEEP_missingHitsThreshold_51"                         , 1   ) ;
  addHEEPParameter(kHC60     , "missingHitsThreshold"                   , "HEEP_missingHitsThreshold_60"                         , 1   ) ;
  
  addHEEPParameter(kHC41     , "dxyFirstPvThresholdBarrel"              , "HEEP_dxyFirstPvThresholdBarrel_41"                    , 0.02) ;
  addHEEPParameter(kHC41     , "dxyFirstPvThresholdEndcap"              , "HEEP_dxyFirstPvThresholdEndcap_41"                    , 0.05) ;
  addHEEPParameter(kHC50_50ns, "dxyFirstPvThresholdBarrel"              , "HEEP_dxyFirstPvThresholdBarrel_50_50ns"               , 0.02) ;
  addHEEPParameter(kHC50_50ns, "dxyFirstPvThresholdEndcap"              , "HEEP_dxyFirstPvThresholdEndcap_50_50ns"               , 0.05) ;
  addHEEPParameter(kHC50_25ns, "dxyFirstPvThresholdBarrel"              , "HEEP_dxyFirstPvThresholdBarrel_50_25ns"               , 0.02) ;
  addHEEPParameter(kHC50_25ns, "dxyFirstPvThresholdEndcap"              , "HEEP_dxyFirstPvThresholdEndcap_50_25ns"               , 0.05) ;
  addHEEPParameter(kHC50     , "dxyFirstPvThresholdBarrel"              , "HEEP_dxyFirstPvThresholdBarrel_50"                    , 0.02) ;
  addHEEPParameter(kHC50     , "dxyFirstPvThresholdEndcap"              , "HEEP_dxyFirstPvThresholdEndcap_50"                    , 0.05) ;
  addHEEPParameter(kHC51     , "dxyFirstPvThresholdBarrel"              , "HEEP_dxyFirstPvThresholdBarrel_51"                    , 0.02) ;
  addHEEPParameter(kHC51     , "dxyFirstPvThresholdEndcap"              , "HEEP_dxyFirstPvThresholdEndcap_51"                    , 0.05) ;
  addHEEPParameter(kHC60     , "dxyFirstPvThresholdBarrel"              , "HEEP_dxyFirstPvThresholdBarrel_60"                    , 0.02) ;
  addHEEPParameter(kHC60     , "dxyFirstPvThresholdEndcap"              , "HEEP_dxyFirstPvThresholdEndcap_60"                    , 0.05) ;
  
  // Isolation
  addHEEPParameter(kHC41     , "isolEMHadDepth1ConstantTermBarrel"      , "HEEP_isolEMHadDepth1ConstantTermBarrel_41"            , 2.0 ) ;
  addHEEPParameter(kHC41     , "isolEMHadDepth1ConstantTermEndcapLowEt" , "HEEP_isolEMHadDepth1ConstantTermEndcapLowEt_41"       , 2.5 ) ;
  addHEEPParameter(kHC41     , "isolEMHadDepth1ConstantTermEndcapHighEt", "HEEP_isolEMHadDepth1ConstantTermEndcapHighEt_41"      , 2.5 ) ;
  addHEEPParameter(kHC41     , "isolEMHadDepth1LinearTermBarrel"        , "HEEP_isolEMHadDepth1LinearTermBarrel_41"              , 0.03) ;
  addHEEPParameter(kHC41     , "isolEMHadDepth1LinearTermEndcap"        , "HEEP_isolEMHadDepth1LinearTermEndcap_41"              , 0.03) ;
  addHEEPParameter(kHC41     , "isolEMHadDepth1OffsetTermEndcap"        , "HEEP_isolEMHadDepth1OffsetTermEndcap_41"              , 50.0) ;
  addHEEPParameter(kHC41     , "isolEMHadDepth1EcalHcal1EffAreaBarrel"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaBarrel_41"        , 0.28) ;
  addHEEPParameter(kHC41     , "isolEMHadDepth1EcalHcal1EffAreaEndcap"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaEndcap_41"        , 0.28) ;
  
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1ConstantTermBarrel"      , "HEEP_isolEMHadDepth1ConstantTermBarrel_50_50ns"       , 2.0 ) ;
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1ConstantTermEndcapLowEt" , "HEEP_isolEMHadDepth1ConstantTermEndcapLowEt_50_50ns"  , 2.5 ) ;
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1ConstantTermEndcapHighEt", "HEEP_isolEMHadDepth1ConstantTermEndcapHighEt_50_50ns" , 2.5 ) ;
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1LinearTermBarrel"        , "HEEP_isolEMHadDepth1LinearTermBarrel_50_50ns"         , 0.03) ;
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1LinearTermEndcap"        , "HEEP_isolEMHadDepth1LinearTermEndcap_50_50ns"         , 0.03) ;
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1OffsetTermEndcap"        , "HEEP_isolEMHadDepth1OffsetTermEndcap_50_50ns"         , 50.0) ;
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1EcalHcal1EffAreaBarrel"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaBarrel_50_50ns"   , 0.28) ;
  addHEEPParameter(kHC50_50ns, "isolEMHadDepth1EcalHcal1EffAreaEndcap"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaEndcap_50_50ns"   , 0.28) ;
  
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1ConstantTermBarrel"      , "HEEP_isolEMHadDepth1ConstantTermBarrel_50_25ns"       , 2.0 ) ;
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1ConstantTermEndcapLowEt" , "HEEP_isolEMHadDepth1ConstantTermEndcapLowEt_50_25ns"  , 2.5 ) ;
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1ConstantTermEndcapHighEt", "HEEP_isolEMHadDepth1ConstantTermEndcapHighEt_50_25ns" , 2.5 ) ;
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1LinearTermBarrel"        , "HEEP_isolEMHadDepth1LinearTermBarrel_50_25ns"         , 0.03) ;
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1LinearTermEndcap"        , "HEEP_isolEMHadDepth1LinearTermEndcap_50_25ns"         , 0.03) ;
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1OffsetTermEndcap"        , "HEEP_isolEMHadDepth1OffsetTermEndcap_50_25ns"         , 50.0) ;
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1EcalHcal1EffAreaBarrel"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaBarrel_50_25ns"   , 0.28) ;
  addHEEPParameter(kHC50_25ns, "isolEMHadDepth1EcalHcal1EffAreaEndcap"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaEndcap_50_25ns"   , 0.28) ;
  
  addHEEPParameter(kHC50     , "isolEMHadDepth1ConstantTermBarrel"      , "HEEP_isolEMHadDepth1ConstantTermBarrel_50"            , 2.0 ) ;
  addHEEPParameter(kHC50     , "isolEMHadDepth1ConstantTermEndcapLowEt" , "HEEP_isolEMHadDepth1ConstantTermEndcapLowEt_50"       , 2.5 ) ;
  addHEEPParameter(kHC50     , "isolEMHadDepth1ConstantTermEndcapHighEt", "HEEP_isolEMHadDepth1ConstantTermEndcapHighEt_50"      , 2.5 ) ;
  addHEEPParameter(kHC50     , "isolEMHadDepth1LinearTermBarrel"        , "HEEP_isolEMHadDepth1LinearTermBarrel_50"              , 0.03) ;
  addHEEPParameter(kHC50     , "isolEMHadDepth1LinearTermEndcap"        , "HEEP_isolEMHadDepth1LinearTermEndcap_50"              , 0.03) ;
  addHEEPParameter(kHC50     , "isolEMHadDepth1OffsetTermEndcap"        , "HEEP_isolEMHadDepth1OffsetTermEndcap_50"              , 50.0) ;
  addHEEPParameter(kHC50     , "isolEMHadDepth1EcalHcal1EffAreaBarrel"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaBarrel_50"        , 0.28) ;
  addHEEPParameter(kHC50     , "isolEMHadDepth1EcalHcal1EffAreaEndcap"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaEndcap_50"        , 0.28) ;
  
  addHEEPParameter(kHC51     , "isolEMHadDepth1ConstantTermBarrel"      , "HEEP_isolEMHadDepth1ConstantTermBarrel_51"            , 2.0 ) ;
  addHEEPParameter(kHC51     , "isolEMHadDepth1ConstantTermEndcapLowEt" , "HEEP_isolEMHadDepth1ConstantTermEndcapLowEt_51"       , 2.5 ) ;
  addHEEPParameter(kHC51     , "isolEMHadDepth1ConstantTermEndcapHighEt", "HEEP_isolEMHadDepth1ConstantTermEndcapHighEt_51"      , 2.5 ) ;
  addHEEPParameter(kHC51     , "isolEMHadDepth1LinearTermBarrel"        , "HEEP_isolEMHadDepth1LinearTermBarrel_51"              , 0.03) ;
  addHEEPParameter(kHC51     , "isolEMHadDepth1LinearTermEndcap"        , "HEEP_isolEMHadDepth1LinearTermEndcap_51"              , 0.03) ;
  addHEEPParameter(kHC51     , "isolEMHadDepth1OffsetTermEndcap"        , "HEEP_isolEMHadDepth1OffsetTermEndcap_51"              , 50.0) ;
  addHEEPParameter(kHC51     , "isolEMHadDepth1EcalHcal1EffAreaBarrel"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaBarrel_51"        , 0.28) ;
  addHEEPParameter(kHC51     , "isolEMHadDepth1EcalHcal1EffAreaEndcap"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaEndcap_51"        , 0.28) ;
  
  addHEEPParameter(kHC60     , "isolEMHadDepth1ConstantTermBarrel"      , "HEEP_isolEMHadDepth1ConstantTermBarrel_60"            , 2.0 ) ;
  addHEEPParameter(kHC60     , "isolEMHadDepth1ConstantTermEndcapLowEt" , "HEEP_isolEMHadDepth1ConstantTermEndcapLowEt_60"       , 2.5 ) ;
  addHEEPParameter(kHC60     , "isolEMHadDepth1ConstantTermEndcapHighEt", "HEEP_isolEMHadDepth1ConstantTermEndcapHighEt_60"      , 2.5 ) ;
  addHEEPParameter(kHC60     , "isolEMHadDepth1LinearTermBarrel"        , "HEEP_isolEMHadDepth1LinearTermBarrel_60"              , 0.03) ;
  addHEEPParameter(kHC60     , "isolEMHadDepth1LinearTermEndcap"        , "HEEP_isolEMHadDepth1LinearTermEndcap_60"              , 0.03) ;
  addHEEPParameter(kHC60     , "isolEMHadDepth1OffsetTermEndcap"        , "HEEP_isolEMHadDepth1OffsetTermEndcap_60"              , 50.0) ;
  addHEEPParameter(kHC60     , "isolEMHadDepth1EcalHcal1EffAreaBarrel"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaBarrel_60"        , 0.28) ;
  addHEEPParameter(kHC60     , "isolEMHadDepth1EcalHcal1EffAreaEndcap"  , "HEEP_isolEMHadDepth1EcalHcal1EffAreaEndcap_60"        , 0.28) ;
  
  addHEEPParameter(kHC41     , "IsolPtTrksThresholdBarrel"              , "HEEP_IsolPtTrksThresholdBarrel_41"                    , 5.0 ) ;
  addHEEPParameter(kHC41     , "IsolPtTrksThresholdEndcap"              , "HEEP_IsolPtTrksThresholdEndcap_41"                    , 5.0 ) ;
  addHEEPParameter(kHC50_50ns, "IsolPtTrksThresholdBarrel"              , "HEEP_IsolPtTrksThresholdBarrel_50_50ns"               , 5.0 ) ;
  addHEEPParameter(kHC50_50ns, "IsolPtTrksThresholdEndcap"              , "HEEP_IsolPtTrksThresholdEndcap_50_50ns"               , 5.0 ) ;
  addHEEPParameter(kHC50_25ns, "IsolPtTrksThresholdBarrel"              , "HEEP_IsolPtTrksThresholdBarrel_50_25ns"               , 5.0 ) ;
  addHEEPParameter(kHC50_25ns, "IsolPtTrksThresholdEndcap"              , "HEEP_IsolPtTrksThresholdEndcap_50_25ns"               , 5.0 ) ;
  addHEEPParameter(kHC50     , "IsolPtTrksThresholdBarrel"              , "HEEP_IsolPtTrksThresholdBarrel_50"                    , 5.0 ) ;
  addHEEPParameter(kHC50     , "IsolPtTrksThresholdEndcap"              , "HEEP_IsolPtTrksThresholdEndcap_50"                    , 5.0 ) ;
  addHEEPParameter(kHC51     , "IsolPtTrksThresholdBarrel"              , "HEEP_IsolPtTrksThresholdBarrel_51"                    , 5.0 ) ;
  addHEEPParameter(kHC51     , "IsolPtTrksThresholdEndcap"              , "HEEP_IsolPtTrksThresholdEndcap_51"                    , 5.0 ) ;
  addHEEPParameter(kHC60     , "IsolPtTrksThresholdBarrel"              , "HEEP_IsolPtTrksThresholdBarrel_60"                    , 5.0 ) ;
  addHEEPParameter(kHC60     , "IsolPtTrksThresholdEndcap"              , "HEEP_IsolPtTrksThresholdEndcap_60"                    , 5.0 ) ;
  
  for(unsigned int i=0 ; i<parameters_.size() ; ++i){
    HEEPParameter* par = parameters_.at(i) ;
    par->setValue(iConfig) ;
  }
}
IIHEModuleHEEP::~IIHEModuleHEEP(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleHEEP::beginJob(){
  std::vector<int> MCPdgIdsToSave ;
  // Leptons
  MCPdgIdsToSave.push_back(11) ; // Electron
  MCPdgIdsToSave.push_back(13) ; // Muon
  MCPdgIdsToSave.push_back(15) ; // Tau
  
  // QCD
  MCPdgIdsToSave.push_back(21) ; // gluon
  MCPdgIdsToSave.push_back( 1) ; // d quark
  MCPdgIdsToSave.push_back( 2) ; // u quark
  MCPdgIdsToSave.push_back( 3) ; // s quark
  MCPdgIdsToSave.push_back( 4) ; // c quark
  MCPdgIdsToSave.push_back( 5) ; // b quark
  MCPdgIdsToSave.push_back( 6) ; // t quark
  
  // SM bosons
  MCPdgIdsToSave.push_back(22) ; // Photon
  MCPdgIdsToSave.push_back(23) ; // Z boson
  MCPdgIdsToSave.push_back(24) ; // W boson
  MCPdgIdsToSave.push_back(25) ; // BEH boson
  
  // BSM boson
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
  
  for(unsigned int i=0 ; i<parameters_.size() ; ++i){
    HEEPParameter* par = parameters_.at(i) ;
    std::string branchName = "HEEP_cutflow_" + par->name() + "_" + par->cutflowName() ;
    addValueToMetaTree(branchName, par->value()) ;
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                     HEEP cutflow                                   //
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                          4.1                                       //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_41 = "HEEP_cutflow41" ;
  HEEPCutflow_41_acceptance_ = new HEEPCutCollection(kHC41, prefix_41 + "_acceptance", this, storeHEEP41_) ;
  HEEPCutflow_41_ID_         = new HEEPCutCollection(kHC41, prefix_41 + "_ID"        , this, storeHEEP41_) ;
  HEEPCutflow_41_isolation_  = new HEEPCutCollection(kHC41, prefix_41 + "_isolation" , this, storeHEEP41_) ;
  HEEPCutflow_41_total_      = new HEEPCutCollection(kHC41, prefix_41 + "_total"     , this, storeHEEP41_) ;
  
  // These cuts must be updated so declare them separately
  cut_41_isolEMHadDepth1_ = new HEEPCut_isolEMHadDepth1(kHC41, prefix_41 + "_isolEMHadDepth1", this) ;
  cut_41_dxyFirstPV_      = new HEEPCut_dxyFirstPV     (kHC41, prefix_41 + "_dxyFirstPV"     , this) ;
  
  HEEPCutflow_41_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_Et (kHC41, prefix_41 + "_Et" , this)) ;
  HEEPCutflow_41_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_eta(kHC41, prefix_41 + "_eta", this)) ;
  
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_EcalDriven   (kHC41, prefix_41 + "_EcalDriven"   , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_dEtaIn    (kHC41, prefix_41 + "_dEtaIn"       , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_dPhiIn       (kHC41, prefix_41 + "_dPhiIn"       , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_41_HOverE    (kHC41, prefix_41 + "_HOverE"       , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_SigmaIetaIeta(kHC41, prefix_41 + "_SigmaIetaIeta", this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_E1x5OverE5x5 (kHC41, prefix_41 + "_E1x5OverE5x5" , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_E2x5OverE5x5 (kHC41, prefix_41 + "_E2x5OverE5x5" , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) new HEEPCut_missingHits  (kHC41, prefix_41 + "_missingHits"  , this)) ;
  HEEPCutflow_41_ID_->addCut( (HEEPCutBase*) cut_41_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_41_isolation_->addCut( (HEEPCutBase*) cut_41_isolEMHadDepth1_) ;
  HEEPCutflow_41_isolation_->addCut( (HEEPCutBase*) new HEEPCut_IsolPtTrks(kHC41, prefix_41 + "_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_41_total_->addCutCollection(HEEPCutflow_41_acceptance_) ;
  HEEPCutflow_41_total_->addCutCollection(HEEPCutflow_41_ID_        ) ;
  HEEPCutflow_41_total_->addCutCollection(HEEPCutflow_41_isolation_ ) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                        5.0 50ns                                    //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_50_50ns = "HEEP_cutflow50_50ns" ;
  HEEPCutflow_50_50ns_acceptance_ = new HEEPCutCollection(kHC50_50ns, prefix_50_50ns+"_acceptance", this, storeHEEP50_50_) ;
  HEEPCutflow_50_50ns_ID_         = new HEEPCutCollection(kHC50_50ns, prefix_50_50ns+"_ID"        , this, storeHEEP50_50_) ;
  HEEPCutflow_50_50ns_isolation_  = new HEEPCutCollection(kHC50_50ns, prefix_50_50ns+"_isolation" , this, storeHEEP50_50_) ;
  HEEPCutflow_50_50ns_total_      = new HEEPCutCollection(kHC50_50ns, prefix_50_50ns+"_total"     , this, storeHEEP50_50_) ;
  
  // These cuts must be updated so declare them separately
  cut_50_50ns_isolEMHadDepth1_ = new HEEPCut_isolEMHadDepth1(kHC50_50ns, prefix_50_50ns + "_isolEMHadDepth1", this) ;
  cut_50_50ns_dxyFirstPV_      = new HEEPCut_dxyFirstPV     (kHC50_50ns, prefix_50_50ns + "_dxyFirstPV"     , this) ;
  
  // Define the ID
  HEEPCutflow_50_50ns_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_Et (kHC50_50ns, prefix_50_50ns + "_Et" , this)) ;
  HEEPCutflow_50_50ns_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_eta(kHC50_50ns, prefix_50_50ns + "_eta", this)) ;
  
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_EcalDriven    (kHC50_50ns, prefix_50_50ns + "_EcalDriven"   , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_dEtaIn(kHC50_50ns, prefix_50_50ns + "_dEtaIn"       , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_dPhiIn        (kHC50_50ns, prefix_50_50ns + "_dPhiIn"       , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_50ns_HOverE(kHC50_50ns, prefix_50_50ns + "_HOverE"       , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_SigmaIetaIeta (kHC50_50ns, prefix_50_50ns + "_SigmaIetaIeta", this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_E1x5OverE5x5  (kHC50_50ns, prefix_50_50ns + "_E1x5OverE5x5" , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_E2x5OverE5x5  (kHC50_50ns, prefix_50_50ns + "_E2x5OverE5x5" , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_missingHits   (kHC50_50ns, prefix_50_50ns + "_missingHits"  , this)) ;
  HEEPCutflow_50_50ns_ID_->addCut( (HEEPCutBase*) cut_50_50ns_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_50_50ns_isolation_->addCut( (HEEPCutBase*) cut_50_50ns_isolEMHadDepth1_) ;
  HEEPCutflow_50_50ns_isolation_->addCut( (HEEPCutBase*) new HEEPCut_IsolPtTrks(kHC50_50ns, prefix_50_50ns + "_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_50_50ns_total_->addCutCollection(HEEPCutflow_50_50ns_acceptance_) ;
  HEEPCutflow_50_50ns_total_->addCutCollection(HEEPCutflow_50_50ns_ID_        ) ;
  HEEPCutflow_50_50ns_total_->addCutCollection(HEEPCutflow_50_50ns_isolation_ ) ;
    
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                        5.0 25ns                                    //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_50_25ns = "HEEP_cutflow50_25ns" ;
  HEEPCutflow_50_25ns_acceptance_ = new HEEPCutCollection(kHC50_25ns, prefix_50_25ns+"_acceptance", this, storeHEEP50_25_) ;
  HEEPCutflow_50_25ns_ID_         = new HEEPCutCollection(kHC50_25ns, prefix_50_25ns+"_ID"        , this, storeHEEP50_25_) ;
  HEEPCutflow_50_25ns_isolation_  = new HEEPCutCollection(kHC50_25ns, prefix_50_25ns+"_isolation" , this, storeHEEP50_25_) ;
  HEEPCutflow_50_25ns_total_      = new HEEPCutCollection(kHC50_25ns, prefix_50_25ns+"_total"     , this, storeHEEP50_25_) ;
  
  // These cuts must be updated so declare them separately
  cut_50_25ns_isolEMHadDepth1_ = new HEEPCut_isolEMHadDepth1(kHC50_25ns, prefix_50_25ns + "_isolEMHadDepth1", this) ;
  cut_50_25ns_dxyFirstPV_      = new HEEPCut_dxyFirstPV     (kHC50_25ns, prefix_50_25ns + "_dxyFirstPV"     , this) ;
  
  // Define the ID
  HEEPCutflow_50_25ns_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_Et (kHC50_25ns, prefix_50_25ns + "_Et" , this)) ;
  HEEPCutflow_50_25ns_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_eta(kHC50_25ns, prefix_50_25ns + "_eta", this)) ;
  
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_EcalDriven    (kHC50_25ns, prefix_50_25ns + "_EcalDriven"   , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_dEtaIn(kHC50_25ns, prefix_50_25ns + "_dEtaIn"       , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_dPhiIn        (kHC50_25ns, prefix_50_25ns + "_dPhiIn"       , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_25ns_HOverE(kHC50_25ns, prefix_50_25ns + "_HOverE"       , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_SigmaIetaIeta (kHC50_25ns, prefix_50_25ns + "_SigmaIetaIeta", this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_E1x5OverE5x5  (kHC50_25ns, prefix_50_25ns + "_E1x5OverE5x5" , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_E2x5OverE5x5  (kHC50_25ns, prefix_50_25ns + "_E2x5OverE5x5" , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) new HEEPCut_missingHits   (kHC50_25ns, prefix_50_25ns + "_missingHits"  , this)) ;
  HEEPCutflow_50_25ns_ID_->addCut( (HEEPCutBase*) cut_50_25ns_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_50_25ns_isolation_->addCut( (HEEPCutBase*) cut_50_25ns_isolEMHadDepth1_) ;
  HEEPCutflow_50_25ns_isolation_->addCut( (HEEPCutBase*) new HEEPCut_IsolPtTrks(kHC50_25ns, prefix_50_25ns + "_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_50_25ns_total_->addCutCollection(HEEPCutflow_50_25ns_acceptance_) ;
  HEEPCutflow_50_25ns_total_->addCutCollection(HEEPCutflow_50_25ns_ID_        ) ;
  HEEPCutflow_50_25ns_total_->addCutCollection(HEEPCutflow_50_25ns_isolation_ ) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                           5.0                                      //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_50 = "HEEP_cutflow50" ;
  HEEPCutflow_50_acceptance_ = new HEEPCutCollection(kHC50, prefix_50+"_acceptance", this, storeHEEP50_) ;
  HEEPCutflow_50_ID_         = new HEEPCutCollection(kHC50, prefix_50+"_ID"        , this, storeHEEP50_) ;
  HEEPCutflow_50_isolation_  = new HEEPCutCollection(kHC50, prefix_50+"_isolation" , this, storeHEEP50_) ;
  HEEPCutflow_50_total_      = new HEEPCutCollection(kHC50, prefix_50+"_total"     , this, storeHEEP50_) ;
  
  // These cuts must be updated so declare them separately
  cut_50_isolEMHadDepth1_ = new HEEPCut_isolEMHadDepth1(kHC50, prefix_50 + "_isolEMHadDepth1", this) ;
  cut_50_dxyFirstPV_      = new HEEPCut_dxyFirstPV     (kHC50, prefix_50 + "_dxyFirstPV"     , this) ;
  
  // Define the ID
  HEEPCutflow_50_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_Et (kHC50, prefix_50 + "_Et" , this)) ;
  HEEPCutflow_50_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_eta(kHC50, prefix_50 + "_eta", this)) ;
  
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_EcalDriven   (kHC50, prefix_50 + "_EcalDriven"   , this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_dEtaIn    (kHC50, prefix_50 + "_dEtaIn"       , this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_dPhiIn       (kHC50, prefix_50 + "_dPhiIn"       , this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_50_HOverE    (kHC50, prefix_50 + "_HOverE"       , this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_SigmaIetaIeta(kHC50, prefix_50 + "_SigmaIetaIeta", this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_E1x5OverE5x5 (kHC50, prefix_50 + "_E1x5OverE5x5" , this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_E2x5OverE5x5 (kHC50, prefix_50 + "_E2x5OverE5x5" , this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) new HEEPCut_missingHits  (kHC50, prefix_50 + "_missingHits"  , this)) ;
  HEEPCutflow_50_ID_->addCut( (HEEPCutBase*) cut_50_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_50_isolation_->addCut( (HEEPCutBase*) cut_50_isolEMHadDepth1_) ;
  HEEPCutflow_50_isolation_->addCut( (HEEPCutBase*) new HEEPCut_IsolPtTrks(kHC50, prefix_50 + "_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_50_total_->addCutCollection(HEEPCutflow_50_acceptance_) ;
  HEEPCutflow_50_total_->addCutCollection(HEEPCutflow_50_ID_        ) ;
  HEEPCutflow_50_total_->addCutCollection(HEEPCutflow_50_isolation_ ) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                           5.1                                      //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_51 = "HEEP_cutflow51" ;
  HEEPCutflow_51_acceptance_ = new HEEPCutCollection(kHC51, prefix_51+"_acceptance", this, storeHEEP51_) ;
  HEEPCutflow_51_ID_         = new HEEPCutCollection(kHC51, prefix_51+"_ID"        , this, storeHEEP51_) ;
  HEEPCutflow_51_isolation_  = new HEEPCutCollection(kHC51, prefix_51+"_isolation" , this, storeHEEP51_) ;
  HEEPCutflow_51_total_      = new HEEPCutCollection(kHC51, prefix_51+"_total"     , this, storeHEEP51_) ;
  
  // These cuts must be updated so declare them separately
  cut_51_isolEMHadDepth1_ = new HEEPCut_isolEMHadDepth1(kHC51, prefix_51 + "_isolEMHadDepth1", this) ;
  cut_51_dxyFirstPV_      = new HEEPCut_dxyFirstPV     (kHC51, prefix_51 + "_dxyFirstPV"     , this) ;
  
  // Define the ID
  HEEPCutflow_51_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_Et (kHC51, prefix_51 + "_Et" , this)) ;
  HEEPCutflow_51_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_eta(kHC51, prefix_51 + "_eta", this)) ;
  
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_EcalDriven   (kHC51, prefix_51 + "_EcalDriven"   , this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_51_dEtaIn    (kHC51, prefix_51 + "_dEtaIn"       , this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_dPhiIn       (kHC51, prefix_51 + "_dPhiIn"       , this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_51_HOverE    (kHC51, prefix_51 + "_HOverE"       , this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_SigmaIetaIeta(kHC51, prefix_51 + "_SigmaIetaIeta", this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_E1x5OverE5x5 (kHC51, prefix_51 + "_E1x5OverE5x5" , this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_E2x5OverE5x5 (kHC51, prefix_51 + "_E2x5OverE5x5" , this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) new HEEPCut_missingHits  (kHC51, prefix_51 + "_missingHits"  , this)) ;
  HEEPCutflow_51_ID_->addCut( (HEEPCutBase*) cut_51_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_51_isolation_->addCut( (HEEPCutBase*) cut_51_isolEMHadDepth1_) ;
  HEEPCutflow_51_isolation_->addCut( (HEEPCutBase*) new HEEPCut_IsolPtTrks(kHC51, prefix_51 + "_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_51_total_->addCutCollection(HEEPCutflow_51_acceptance_) ;
  HEEPCutflow_51_total_->addCutCollection(HEEPCutflow_51_ID_        ) ;
  HEEPCutflow_51_total_->addCutCollection(HEEPCutflow_51_isolation_ ) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                                           6.0                                      //
  ////////////////////////////////////////////////////////////////////////////////////////
  std::string prefix_60 = "HEEP_cutflow60" ;
  HEEPCutflow_60_acceptance_ = new HEEPCutCollection(kHC60, prefix_60+"_acceptance", this, storeHEEP50_) ;
  HEEPCutflow_60_ID_         = new HEEPCutCollection(kHC60, prefix_60+"_ID"        , this, storeHEEP50_) ;
  HEEPCutflow_60_isolation_  = new HEEPCutCollection(kHC60, prefix_60+"_isolation" , this, storeHEEP50_) ;
  HEEPCutflow_60_total_      = new HEEPCutCollection(kHC60, prefix_60+"_total"     , this, storeHEEP50_) ;
  
  // These cuts must be updated so declare them separately
  cut_60_isolEMHadDepth1_ = new HEEPCut_isolEMHadDepth1(kHC60, prefix_60 + "_isolEMHadDepth1", this) ;
  cut_60_dxyFirstPV_      = new HEEPCut_dxyFirstPV     (kHC60, prefix_60 + "_dxyFirstPV"     , this) ;
  
  // Define the ID
  HEEPCutflow_60_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_Et (kHC60, prefix_60 + "_Et" , this)) ;
  HEEPCutflow_60_acceptance_->addCut( (HEEPCutBase*) new HEEPCut_eta(kHC60, prefix_60 + "_eta", this)) ;
  
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_EcalDriven   (kHC60, prefix_60 + "_EcalDriven"   , this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_60_dEtaIn    (kHC60, prefix_60 + "_dEtaIn"       , this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_dPhiIn       (kHC60, prefix_60 + "_dPhiIn"       , this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_60_HOverE    (kHC60, prefix_60 + "_HOverE"       , this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_SigmaIetaIeta(kHC60, prefix_60 + "_SigmaIetaIeta", this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_E1x5OverE5x5 (kHC60, prefix_60 + "_E1x5OverE5x5" , this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_E2x5OverE5x5 (kHC60, prefix_60 + "_E2x5OverE5x5" , this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) new HEEPCut_missingHits  (kHC60, prefix_60 + "_missingHits"  , this)) ;
  HEEPCutflow_60_ID_->addCut( (HEEPCutBase*) cut_60_dxyFirstPV_) ;
  
  // Define the isolation
  HEEPCutflow_60_isolation_->addCut( (HEEPCutBase*) cut_60_isolEMHadDepth1_) ;
  HEEPCutflow_60_isolation_->addCut( (HEEPCutBase*) new HEEPCut_IsolPtTrks(kHC60, prefix_60 + "_IsolPtTrks", this)) ;
  
  // Put it all together
  HEEPCutflow_60_total_->addCutCollection(HEEPCutflow_60_acceptance_) ;
  HEEPCutflow_60_total_->addCutCollection(HEEPCutflow_60_ID_        ) ;
  HEEPCutflow_60_total_->addCutCollection(HEEPCutflow_60_isolation_ ) ;
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //                     Add everything to the vector of cutflows                       //
  ////////////////////////////////////////////////////////////////////////////////////////
  HEEPCutflows_.push_back(HEEPCutflow_41_acceptance_     ) ;
  HEEPCutflows_.push_back(HEEPCutflow_41_ID_             ) ;
  HEEPCutflows_.push_back(HEEPCutflow_41_isolation_      ) ;
  HEEPCutflows_.push_back(HEEPCutflow_41_total_          ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_25ns_acceptance_) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_25ns_ID_        ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_25ns_isolation_ ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_25ns_total_     ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_50ns_acceptance_) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_50ns_ID_        ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_50ns_isolation_ ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_50ns_total_     ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_acceptance_     ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_ID_             ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_isolation_      ) ;
  HEEPCutflows_.push_back(HEEPCutflow_50_total_          ) ;
  HEEPCutflows_.push_back(HEEPCutflow_51_acceptance_     ) ;
  HEEPCutflows_.push_back(HEEPCutflow_51_ID_             ) ;
  HEEPCutflows_.push_back(HEEPCutflow_51_isolation_      ) ;
  HEEPCutflows_.push_back(HEEPCutflow_51_total_          ) ;
  HEEPCutflows_.push_back(HEEPCutflow_60_ID_             ) ;
  HEEPCutflows_.push_back(HEEPCutflow_60_isolation_      ) ;
  HEEPCutflows_.push_back(HEEPCutflow_60_total_          ) ;
  
  for(unsigned int i=0 ; i<HEEPCutflows_.size() ; ++i){
    if(HEEPCutflows_.at(i)->isActive()==false) continue ;
    HEEPCutflows_.at(i)->config(parameters_) ;
  }
}

// ------------ method called to for each event  ------------
void IIHEModuleHEEP::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
  
  // Pass parameters to the HEEP cutflow objects
  math::XYZPoint* firstPrimaryVertex = parent_->getFirstPrimaryVertex() ;
  cut_41_dxyFirstPV_     ->setFirstPV(firstPrimaryVertex) ;
  cut_50_50ns_dxyFirstPV_->setFirstPV(firstPrimaryVertex) ;
  cut_50_25ns_dxyFirstPV_->setFirstPV(firstPrimaryVertex) ;
  cut_50_dxyFirstPV_     ->setFirstPV(firstPrimaryVertex) ;
  cut_51_dxyFirstPV_     ->setFirstPV(firstPrimaryVertex) ;
  cut_60_dxyFirstPV_     ->setFirstPV(firstPrimaryVertex) ;
  
  edm::Handle<double> rhoHandle ;
  iEvent.getByLabel(InputTag("fixedGridRhoAll"), rhoHandle) ;
  double rho = *rhoHandle ;
  cut_41_isolEMHadDepth1_     ->setRho(rho) ;
  cut_50_50ns_isolEMHadDepth1_->setRho(rho) ;
  cut_50_25ns_isolEMHadDepth1_->setRho(rho) ;
  cut_50_isolEMHadDepth1_     ->setRho(rho) ;
  cut_51_isolEMHadDepth1_     ->setRho(rho) ;
  cut_60_isolEMHadDepth1_     ->setRho(rho) ;
  
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
  
/*CHOOSE_RELEASE_START DEFAULT
    EcalClusterLazyTools lazytool(iEvent, iSetup, parent_->getReducedBarrelRecHitCollectionToken(), parent_->getReducedEndcapRecHitCollectionToken(), parent_->getReducedESRecHitCollectionToken()) ;
CHOOSE_RELEASE_END DEFAULT*/
//CHOOSE_RELEASE_START CMSSW_7_3_0 CMSSW_7_2_0
    EcalClusterLazyTools lazytool(iEvent, iSetup, parent_->getReducedBarrelRecHitCollectionToken(), parent_->getReducedEndcapRecHitCollectionToken()) ;
//CHOOSE_RELEASE_END CMSSW_7_3_0 CMSSW_7_2_0
/*CHOOSE_RELEASE_START CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11
    EcalClusterLazyTools lazytool(iEvent, iSetup, InputTag("reducedEcalRecHitsEB"),InputTag("reducedEcalRecHitsEE"),InputTag("reducedEcalRecHitsES"));
CHOOSE_RELEASE_END CMSSW_7_0_6_patch1 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11*/
  
  reco::GsfElectronCollection electrons = parent_->getElectronCollection() ;
  for(reco::GsfElectronCollection::const_iterator gsfiter = electrons.begin() ; gsfiter!=electrons.end() ; ++gsfiter){
    reco::GsfElectron* gsf = (reco::GsfElectron*) &*gsfiter ;
    
    if(storeHEEP41_   ) HEEPCutflow_41_total_     ->applyCuts(gsf) ;
    if(storeHEEP50_50_) HEEPCutflow_50_50ns_total_->applyCuts(gsf) ;
    if(storeHEEP50_25_) HEEPCutflow_50_25ns_total_->applyCuts(gsf) ;
    if(storeHEEP50_   ) HEEPCutflow_50_total_     ->applyCuts(gsf) ;
    if(storeHEEP51_   ) HEEPCutflow_51_total_     ->applyCuts(gsf) ;
    if(storeHEEP60_   ) HEEPCutflow_60_total_     ->applyCuts(gsf) ;
    
    // Required for preshower variables
    reco::SuperClusterRef    cl_ref = gsf->superCluster() ;
    const reco::CaloClusterPtr seed = gsf->superCluster()->seed() ;
    
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
  bool accept = false ;
  for(unsigned i=0 ; i<HEEPCutflows_.size() ; i++){
    HEEPCutflows_.at(i)->endEvent() ;
    if(HEEPCutflows_.at(i)->nPass()>0){
      accept = true ;
      break ;
    }
  }
  if(accept){
    acceptEvent() ;
    nAccept_++ ;
  }
}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleHEEP::endJob(){
  std::cout << std::endl << "IIHEModuleHEEP report:" << std::endl ;
  std::cout << "  nAccept  = " << nAccept_   << std::endl ;
}

DEFINE_FWK_MODULE(IIHEModuleHEEP);
