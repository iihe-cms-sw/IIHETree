##########################################################################################
#                                     PHYS14 samples                                     #
##########################################################################################
# Signal
add_job('PHYS14_ZprimeToEE_M5000_20bx25', '/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM')
add_job('PHYS14_ZprimeToEE_M5000_30bx50', '/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/Phys14DR-PU30bx50_PHYS14_25_V1-v1/AODSIM')
add_job('PHYS14_ZprimeToEE_M5000_40bx25', '/ZprimeToEE_M-5000_Tune4C_13TeV-pythia8/Phys14DR-PU40bx25_PHYS14_25_V1-v1/AODSIM')

# Drell-Yan
add_job('PHYS14_DYToEE_20BX25'        , '/DYToEE_Tune4C_13TeV-pythia8/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_DYToEE_30BX50'        , '/DYToEE_Tune4C_13TeV-pythia8/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_DYToEE_Flat20To50BX50', '/DYToEE_Tune4C_13TeV-pythia8/Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_DYToMM_20BX25'        , '/DYToMuMu_Tune4C_13TeV-pythia8/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_DYToMM_30BX50'        , '/DYToMuMu_Tune4C_13TeV-pythia8/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_DYToMM_Flat20To50BX50', '/DYToMuMu_Tune4C_13TeV-pythia8/Phys14DR-Flat20to50BX50_tsg_PHYS14_ST_V1-v1/AODSIM')

# Electroweak (incomplete)
add_job('PHYS14_TT_AVE30BX50'         , '/TT_Tune4C_13TeV-pythia8-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_TT_20bx25'            , '/TT_Tune4C_13TeV-pythia8-tauola/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM')
add_job('PHYS14_WEnu_20bx25'          , 'WToENu_Tune4C_13TeV-pythia8/Phys14DR-PU20bx25_tsg_PHYS14_25_V1-v1/AODSIM')

# QCD (incomplete)
add_job('PHYS14_QCD_50_80_20bx25'    , '/QCD_Pt-50to80_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM'  )
add_job('PHYS14_QCD_50_80_30bx50'    , '/QCD_Pt-50to80_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM'  )
add_job('PHYS14_QCD_80_120_20bx25'   , '/QCD_Pt-80to120_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM' )
add_job('PHYS14_QCD_80_120_30bx50'   , '/QCD_Pt-80to120_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM' )
add_job('PHYS14_QCD_120_170_20bx25'  , '/QCD_Pt-120to170_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_QCD_120_170_30bx50'  , '/QCD_Pt-120to170_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM')

add_job('PHYS14_QCD_470_600_20bx25'  , '/QCD_Pt-470to600_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_QCD_470_600_30bx50'  , '/QCD_Pt-470to600_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_QCD_600_800_20bx25'  , '/QCD_Pt-600to800_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_QCD_600_800_30bx50'  , '/QCD_Pt-600to800_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_QCD_800_1000_20bx25' , '/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_QCD_800_1000_30bx50' , '/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_QCD_1000_1400_20bx25', '/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_QCD_1000_1400_30bx50', '/QCD_Pt-1000to1400_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_QCD_1400_1800_20bx25', '/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_QCD_1400_1800_30bx50', '/QCD_Pt-1400to1800_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM')
add_job('PHYS14_QCD_1800_20bx25'     , '/QCD_Pt-1800_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/AODSIM')
add_job('PHYS14_QCD_1800_30bx50'     , '/QCD_Pt-1800_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/AODSIM')
