#ifndef UserCode_IIHETree_Types_h
#define UserCode_IIHETree_Types_h

enum variableTypes{
  kBool,
  kDouble,
  kFloat,
  kInt,
  kUInt,
  kVectorBool,
  kVectorDouble,
  kVectorFloat,
  kVectorInt,
  kVectorUInt,
  kVectorVectorBool,
  kVectorVectorDouble,
  kVectorVectorFloat,
  kVectorVectorInt,
  kVectorVectorUInt
};

enum particleTypes{
  kSuperCluster,
  kPhoton,
  kElectron,
  kMuon,
  kTau,
  kJet,
  kMET
};

enum triggerLevels{
  kLevel1,
  kHighLevel
};

#endif
