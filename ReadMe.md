Readme


This is currently under construction, 

# Setting up the environment:

```bash
$ setenv SCRAM_ARCH el8_amd64_gcc10 

$ cmsrel CMSSW_12_6_0

$ cd CMSSW_12_6_0/src

$ cmsenv

$ git cms-init

$ git clone --recursive  https://github.com/steffi-bower/PNet_Skimmer.git

$ git cms-addpkg PhysicsTools/PatAlgos

$ rm PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc ##Replace with my own writen by redwan

$ mv PNetSkimmer/PATTauSlimmer.cc PhysicsTools/PatAlgos/plugins/

$ scram b -j 8
```


# Running the code:

```bash
$ cd CMSSW_12_6_0/src

$ scram b clean

$ scram b -j 8

$ cd test/

$ cmsRun rerunTauRecoOnMiniAOD_WithClean_Custom.py
```
