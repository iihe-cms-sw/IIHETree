#ifndef UserCode_IIHETree_IIHEAnalysis_h
#define UserCode_IIHETree_IIHEAnalysis_h

// system include files
#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "UserCode/IIHETree/interface/BranchWrapper.h"
#include "UserCode/IIHETree/interface/IIHEModule.h"

#include "TFile.h"
#include "TTree.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

enum variableTypes{
  kBool,
  kDouble,
  kFloat,
  kInt,
  kVectorBool,
  kVectorDouble,
  kVectorFloat,
  kVectorInt,
  kVectorVectorBool,
  kVectorVectorDouble,
  kVectorVectorFloat,
  kVectorVectorInt
};


class IIHEModule ; // Forward declaration

// class decleration
class IIHEAnalysis : public edm::EDAnalyzer {

private:

public:
  explicit IIHEAnalysis(const edm::ParameterSet& iConfig);
  ~IIHEAnalysis();
  
  bool store(std::string, bool    );
  bool store(std::string, double  );
  bool store(std::string, float   );
  bool store(std::string, int     );
  bool store(std::string, unsigned);
  bool store(std::string, std::vector<bool>  );
  bool store(std::string, std::vector<double>);
  bool store(std::string, std::vector<float >);
  bool store(std::string, std::vector<int   >);
  
  bool addBranch(std::string);
  bool addBranch(std::string,int);
  bool branchExists(std::string);
  
  void setBranchType(int);
  int  getBranchType();
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  
  void beginEvent() ;
  void endEvent() ;
  
  void configureBranches();
  
  // ----------member data ---------------------------
  std::vector<BranchWrapperBase*> allVars_ ;
  std::vector<BranchWrapperBVV* > vars_BVV_;
  std::vector<BranchWrapperDVV* > vars_DVV_;
  std::vector<BranchWrapperFVV* > vars_FVV_;
  std::vector<BranchWrapperIVV* > vars_IVV_;
  std::vector<BranchWrapperBV*  > vars_BV_ ;
  std::vector<BranchWrapperDV*  > vars_DV_ ;
  std::vector<BranchWrapperFV*  > vars_FV_ ;
  std::vector<BranchWrapperIV*  > vars_IV_ ;
  std::vector<BranchWrapperB*   > vars_B_  ;
  std::vector<BranchWrapperD*   > vars_D_  ;
  std::vector<BranchWrapperF*   > vars_F_  ;
  std::vector<BranchWrapperI*   > vars_I_  ;
  
  int currentVarType_ ;
  std::vector< std::pair<std::string, int> > listOfBranches_ ;
  std::vector< std::pair<std::string, int> > missingBranches_ ;
  
  bool debug_;

  // config parameters -------------------------------
  TFile* myFile ;
  TTree* mytree ;
  
  std::vector<IIHEModule*> childModules_;
};
#endif
//define this as a plug-in

