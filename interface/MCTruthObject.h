#ifndef UserCode_IIHETree_MCTruthObject_h
#define UserCode_IIHETree_MCTruthObject_h

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Candidate/interface/Candidate.h"

class MCTruthObject{
private:
  const reco::Candidate* candidate_ ;
  float DeltaRCut_ ;
  std::vector<const reco::Candidate*> mothers_ ;
public:
  MCTruthObject(reco::Candidate*) ;
  ~MCTruthObject() ;
  void addMother(const reco::Candidate*) ;
  int matchMother(std::vector<MCTruthObject*>, unsigned int) ;
  const reco::Candidate* getCandidate(){ return candidate_ ; }
  const reco::Candidate* getMother(unsigned int) ;
  unsigned nMothers(){ return mothers_.size() ; }
};

#endif
