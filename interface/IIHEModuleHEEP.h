#ifndef UserCode_IIHETree_IIHEModuleHEEP_h
#define UserCode_IIHETree_IIHEModuleHEEP_h

class HEEPCut_isolEMHadDepth1 ;
class HEEPCut_dxyFirstPV ;
class HEEPCutCollection ;

#include "UserCode/IIHETree/interface/IIHEModule.h"
#include "UserCode/IIHETree/interface/HEEPCut.h"

// class decleration
class IIHEModuleHEEP : public IIHEModule {
private:
  bool storeHEEP41_    ;
  bool storeHEEP50_50_ ;
  bool storeHEEP50_25_ ;
  bool storeHEEP50_    ;
  bool storeHEEP51_    ;

  double EcalHcal1EffAreaBarrel_  ;
  double EcalHcal1EffAreaEndcaps_ ;
  double rho_ ;
  
  // Define the ID
  
  // Acceptance
  float EtThresholdBarrel_41_ ;
  float EtThresholdEndcap_41_ ;
  float EtThresholdBarrel_50_ ;
  float EtThresholdEndcap_50_ ;
  
  float barrelEtaUpper_41_ ;
  float endcapEtaLower_41_ ;
  float endcapEtaUpper_41_ ;
  float barrelEtaUpper_50_ ;
  float endcapEtaLower_50_ ;
  float endcapEtaUpper_50_ ;
  
  // ID
  float dEtaInThresholdBarrel_41_ ;
  float dEtaInThresholdEndcap_41_ ;
  
  float dEtaInConstantTermBarrel_50_50ns_ ;
  float dEtaInLinearTermBarrel_50_50ns_   ;
  float dEtaInCutoffTermBarrel_50_50ns_   ;
  float dEtaInThresholdEndcap_50_50ns_    ;
  
  float dEtaInConstantTermBarrel_50_25ns_ ;
  float dEtaInLinearTermBarrel_50_25ns_   ;
  float dEtaInCutoffTermBarrel_50_25ns_   ;
  float dEtaInConstantTermEndcap_50_25ns_ ;
  float dEtaInLinearTermEndcap_50_25ns_   ;
  float dEtaInCutoffTermEndcap_50_25ns_   ;
  
  float dEtaInConstantTermBarrel_50_ ;
  float dEtaInLinearTermBarrel_50_   ;
  float dEtaInCutoffTermBarrel_50_   ;
  float dEtaInConstantTermEndcap_50_ ;
  float dEtaInLinearTermEndcap_50_   ;
  float dEtaInCutoffTermEndcap_50_   ;
  
  float dEtaInThresholdBarrel_51_ ;
  float dEtaInThresholdEndcap_51_ ;
  
  float dPhiInThresholdBarrel_41_      ;
  float dPhiInThresholdEndcap_41_      ;
  float dPhiInThresholdBarrel_50_50ns_ ;
  float dPhiInThresholdEndcap_50_50ns_ ;
  float dPhiInThresholdBarrel_50_25ns_ ;
  float dPhiInThresholdEndcap_50_25ns_ ;
  float dPhiInThresholdBarrel_50_      ;
  float dPhiInThresholdEndcap_50_      ;
  float dPhiInThresholdBarrel_51_      ;
  float dPhiInThresholdEndcap_51_      ;
  
  float HOverEThresholdBarrel_41_ ;
  float HOverEThresholdEndcap_41_ ;
  
  float HOverEReciprocalTermBarrel_50_50ns_ ;
  float HOverEConstantTermBarrel_50_50ns_   ;
  float HOverEReciprocalTermEndcap_50_50ns_ ;
  float HOverEConstantTermEndcap_50_50ns_   ;
  
  float HOverEReciprocalTermBarrel_50_25ns_ ;
  float HOverEConstantTermBarrel_50_25ns_   ;
  float HOverEReciprocalTermEndcap_50_25ns_ ;
  float HOverEConstantTermEndcap_50_25ns_   ;
  
  float HOverEReciprocalTermBarrel_50_ ;
  float HOverEConstantTermBarrel_50_   ;
  float HOverEReciprocalTermEndcap_50_ ;
  float HOverEConstantTermEndcap_50_   ;
  
  float HOverEReciprocalTermBarrel_51_ ;
  float HOverEConstantTermBarrel_51_   ;
  float HOverEReciprocalTermEndcap_51_ ;
  float HOverEConstantTermEndcap_51_   ;
  
  float SigmaIetaIetaThreshold_41_      ;
  float SigmaIetaIetaThreshold_50_50ns_ ;
  float SigmaIetaIetaThreshold_50_25ns_ ;
  float SigmaIetaIetaThreshold_50_      ;
  float SigmaIetaIetaThreshold_51_      ;
  
  float E1x5threshold_41_      ;
  float E2x5threshold_41_      ;
  float E1x5threshold_50_50ns_ ;
  float E2x5threshold_50_50ns_ ;
  float E1x5threshold_50_25ns_ ;
  float E2x5threshold_50_25ns_ ;
  float E1x5threshold_50_      ;
  float E2x5threshold_50_      ;
  float E1x5threshold_51_      ;
  float E2x5threshold_51_      ;
  
  float missingHitsThreshold_41_      ;
  float missingHitsThreshold_50_50ns_ ;
  float missingHitsThreshold_50_25ns_ ;
  float missingHitsThreshold_50_      ;
  float missingHitsThreshold_51_      ;
  
  float dxyFirstPvThresholdBarrel_41_      ;
  float dxyFirstPvThresholdEndcap_41_      ;
  float dxyFirstPvThresholdBarrel_50_50ns_ ;
  float dxyFirstPvThresholdEndcap_50_50ns_ ;
  float dxyFirstPvThresholdBarrel_50_25ns_ ;
  float dxyFirstPvThresholdEndcap_50_25ns_ ;
  float dxyFirstPvThresholdBarrel_50_      ;
  float dxyFirstPvThresholdEndcap_50_      ;
  float dxyFirstPvThresholdBarrel_51_      ;
  float dxyFirstPvThresholdEndcap_51_      ;
  
  // Isolation
  float isolEMHadDepth1ConstantTermBarrel_41_       ;
  float isolEMHadDepth1ConstantTermEndcapLowEt_41_  ;
  float isolEMHadDepth1ConstantTermEndcapHighEt_41_ ;
  float isolEMHadDepth1LinearTermBarrel_41_         ;
  float isolEMHadDepth1LinearTermEndcap_41_         ;
  float isolEMHadDepth1OffsetTermEndcap_41_         ;
  
  float isolEMHadDepth1ConstantTermBarrel_50_50ns_       ;
  float isolEMHadDepth1ConstantTermEndcapLowEt_50_50ns_  ;
  float isolEMHadDepth1ConstantTermEndcapHighEt_50_50ns_ ;
  float isolEMHadDepth1LinearTermBarrel_50_50ns_         ;
  float isolEMHadDepth1LinearTermEndcap_50_50ns_         ;
  float isolEMHadDepth1OffsetTermEndcap_50_50ns_         ;
  
  float isolEMHadDepth1ConstantTermBarrel_50_25ns_       ;
  float isolEMHadDepth1ConstantTermEndcapLowEt_50_25ns_  ;
  float isolEMHadDepth1ConstantTermEndcapHighEt_50_25ns_ ;
  float isolEMHadDepth1LinearTermBarrel_50_25ns_         ;
  float isolEMHadDepth1LinearTermEndcap_50_25ns_         ;
  float isolEMHadDepth1OffsetTermEndcap_50_25ns_         ;
  
  float isolEMHadDepth1ConstantTermBarrel_50_       ;
  float isolEMHadDepth1ConstantTermEndcapLowEt_50_  ;
  float isolEMHadDepth1ConstantTermEndcapHighEt_50_ ;
  float isolEMHadDepth1LinearTermBarrel_50_         ;
  float isolEMHadDepth1LinearTermEndcap_50_         ;
  float isolEMHadDepth1OffsetTermEndcap_50_         ;
  
  float isolEMHadDepth1ConstantTermBarrel_51_       ;
  float isolEMHadDepth1ConstantTermEndcapLowEt_51_  ;
  float isolEMHadDepth1ConstantTermEndcapHighEt_51_ ;
  float isolEMHadDepth1LinearTermBarrel_51_         ;
  float isolEMHadDepth1LinearTermEndcap_51_         ;
  float isolEMHadDepth1OffsetTermEndcap_51_         ;
  
  float IsolPtTrksThresholdBarrel_41_      ;
  float IsolPtTrksThresholdEndcap_41_      ;
  float IsolPtTrksThresholdBarrel_50_50ns_ ;
  float IsolPtTrksThresholdEndcap_50_50ns_ ;
  float IsolPtTrksThresholdBarrel_50_25ns_ ;
  float IsolPtTrksThresholdEndcap_50_25ns_ ;
  float IsolPtTrksThresholdBarrel_50_      ;
  float IsolPtTrksThresholdEndcap_50_      ;
  float IsolPtTrksThresholdBarrel_51_      ;
  float IsolPtTrksThresholdEndcap_51_      ;
  
  HEEPCut_isolEMHadDepth1* cut_41_isolEMHadDepth1_ ;
  HEEPCut_dxyFirstPV*      cut_41_dxyFirstPV_      ;
  HEEPCutCollection* HEEPCutflow_41_acceptance_ ;
  HEEPCutCollection* HEEPCutflow_41_ID_         ;
  HEEPCutCollection* HEEPCutflow_41_isolation_  ;
  HEEPCutCollection* HEEPCutflow_41_total_      ;
  
  HEEPCut_isolEMHadDepth1* cut_50_50ns_isolEMHadDepth1_ ;
  HEEPCut_dxyFirstPV*      cut_50_50ns_dxyFirstPV_      ;
  HEEPCutCollection* HEEPCutflow_50_50ns_acceptance_ ;
  HEEPCutCollection* HEEPCutflow_50_50ns_ID_         ;
  HEEPCutCollection* HEEPCutflow_50_50ns_isolation_  ;
  HEEPCutCollection* HEEPCutflow_50_50ns_total_      ;
  
  HEEPCut_isolEMHadDepth1* cut_50_25ns_isolEMHadDepth1_ ;
  HEEPCut_dxyFirstPV*      cut_50_25ns_dxyFirstPV_      ;
  HEEPCutCollection* HEEPCutflow_50_25ns_acceptance_ ;
  HEEPCutCollection* HEEPCutflow_50_25ns_ID_         ;
  HEEPCutCollection* HEEPCutflow_50_25ns_isolation_  ;
  HEEPCutCollection* HEEPCutflow_50_25ns_total_      ;
  
  HEEPCut_isolEMHadDepth1* cut_50_isolEMHadDepth1_ ;
  HEEPCut_dxyFirstPV*      cut_50_dxyFirstPV_      ;
  HEEPCutCollection* HEEPCutflow_50_acceptance_ ;
  HEEPCutCollection* HEEPCutflow_50_ID_         ;
  HEEPCutCollection* HEEPCutflow_50_isolation_  ;
  HEEPCutCollection* HEEPCutflow_50_total_      ;
  
  HEEPCut_isolEMHadDepth1* cut_51_isolEMHadDepth1_ ;
  HEEPCut_dxyFirstPV*      cut_51_dxyFirstPV_      ;
  HEEPCutCollection* HEEPCutflow_51_acceptance_ ;
  HEEPCutCollection* HEEPCutflow_51_ID_         ;
  HEEPCutCollection* HEEPCutflow_51_isolation_  ;
  HEEPCutCollection* HEEPCutflow_51_total_      ;
  
  std::vector<std::string> triggersForMatching_ ;
public:
  explicit IIHEModuleHEEP(const edm::ParameterSet& iConfig);
  ~IIHEModuleHEEP();
  
  void   pubBeginJob(){   beginJob() ; } ;
  void pubBeginEvent(){ beginEvent() ; } ;
  void   pubEndEvent(){   endEvent() ; } ;
  virtual void pubAnalyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){ analyze(iEvent, iSetup) ; } ;
  
  virtual void beginEvent() ;
  virtual void endEvent() ;
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  
  std::vector<HEEPCutCollection*> HEEPCutflows_ ;
};

#endif
