#!/bin/bash

source $VO_CMS_SW_DIR/cmsset_default.sh

dir=`pwd`
if [[ ${dir} == *"CMSSW"* ]]
then
  echo "You are in a CMSSW release.  This script normally checks out different CMSSW releases, so to avoid conflicts you cannot execute it from inside a CMSSW release."
  echo "I will quit now.  Have a good day."
  exit 0
fi

workingDir=IIHETree_testBuild
mkdir ${workingDir}
cd ${workingDir}
for release in CMSSW_7_0_6_patch1 CMSSW_7_3_0 CMSSW_7_2_0 CMSSW_6_2_5 CMSSW_6_2_0_SLHC23_patch1 CMSSW_5_3_11
do
    if [[ ${release} == "CMSSW_7_0_6_patch1" ]]
    then
      export SCRAM_ARCH=slc6_amd64_gcc481
    fi
    
    if [[ ${release} == "CMSSW_7_3_0" ]]
    then
      export SCRAM_ARCH=slc6_amd64_gcc491
    fi
    
    if [[ ${release} == "CMSSW_6_2_5" ]]
    then
      export SCRAM_ARCH=slc5_amd64_gcc472
    fi
    
    if [[ ${release} == "CMSSW_6_2_0_SLHC23_patch1" ]]
    then
      export SCRAM_ARCH=slc6_amd64_gcc472
    fi
    
    if [[ ${release} == "CMSSW_5_3_11" ]]
    then
      export SCRAM_ARCH=slc5_amd64_gcc462
    fi
    
    echo "Working with release "${release}
    scram project CMSSW ${release}
    cd ${release}/src
    eval `scramv1 runtime -sh`
    cp -r ../../../UserCode .
    cp UserCode/IIHETree/test/change_release.py .
    python change_release.py --release=${release}
    scram b clean
    scram b -j4
    echo
    echo
    cd ../../
done

