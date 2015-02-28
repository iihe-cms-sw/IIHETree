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

// These are used in triggers
enum particleTypes{
  kSuperCluster,
  kPhoton,
  kElectron,
  kMuon,
  kTau,
  kJet,
  kBJet,
  kMET,
  kHT,
  kALCa
};

enum triggerLevels{
  kLevel1,
  kHighLevel
};

#endif
