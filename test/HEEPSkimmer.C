#include <TString.h>
#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <vector>

using namespace std ;

void HEEPSkimmer(TString input_name="outfile.root", TString output_name="outfile_skimmed.root", TString cut="1==1"){
  TFile* file_in = new TFile(input_name,"READ") ;
  TTree* tree_in = (TTree*)file_in->Get("IIHEAnalysis") ;
  
  // Select only the branches we want to save
  vector<TString> branchPatterns ;
  branchPatterns.push_back("ev_*"  ) ;
  branchPatterns.push_back("sc_*"  ) ;
  branchPatterns.push_back("gsf_*" ) ;
  branchPatterns.push_back("HEEP_*") ;
  
  tree_in->SetBranchStatus("*",0) ;
  for(unsigned int i=0 ; i<branchPatterns.size() ; i++){
    tree_in->SetBranchStatus(branchPatterns.at(i) ,1) ;
  }
  
  // Finally, copy the tree with some simple cuts
  TFile* file_out = new TFile(output_name, "RECREATE") ;
  TTree* tree_out = tree_in->CopyTree(cut) ;
  
  file_out->Write() ;
}


