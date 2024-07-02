# ECAL CLuatering for Boosted diphotons using CLUE

Installation and complilation:
- Requires `cmssw-el7`
```
cmsrel CMSSW_13_2_4
cd CMSSW_13_2_4/src/
cmssw-el7 -- cmsenv
mkdir test
cd test
git clone https://gitlab.cern.ch/gagarwal/clusteringanalyzer.git
cd ClusteringAnalyzer
cmssw-el7 -- scram b
```
This script is made to take a list of ECAL crystals and PS deposits from MiniAOD and cluster them using CLUE algorithm.
The main script being called by the python config is `plugins/ClusteringAnalyzer.cc`.
To run:
```
cmssw-el7
cmsenv
cmsRun python/ConfFile_cfg.py >& log.txt & # run and safe the stdout to log.txt
```
