Readme


This is currently under construction, 

# Setting up the environment:

```bash
$ setenv SCRAM_ARCH el8_amd64_gcc10 

$ cmsrel CMSSW_12_6_0

$ cd CMSSW_12_6_0/src

$ git cms-init

$ cmsenv

$ git clone --recursive git@github.com:red1habibullah/MiniAODSkimmer.git -b UL_12X_2018

$ git cms-addpkg PhysicsTools/PatAlgos

$ rm PhysicsTools/PatAlgos/plugins/PATTauSlimmer.cc ##Replace with my own writen by redwan


$ scram b -j 8
```


# Running the code:

```bash
$ cd CMSSW_12_6_0/src

$ scram b clean

$ scram b -j8

$ cd test/

$ cmsRun rerunTauRecoOnMiniAOD_WithClean_Custom.py
```
