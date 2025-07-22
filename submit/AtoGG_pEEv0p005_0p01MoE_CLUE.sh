#!/bin/bash

inputfile=/eos/user/g/gziemyte/AtoGG_500events0p005_0p01MoE10p0_200p0EPileupNoResv1_MiniAODpEE.root
version=pEE
MoE=0p005_0p01
active_directory=/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/EEdc1rhoc15delt20
file_name="hist_AtoGG_5000events_${MoE}_Ma2dc15rhoc3deltnoresv${version}.root"

echo "CLUE Clustering"
cd $active_directory
eval `scramv1 runtime -sh`
cmsRun ../python/ConfFile_cfg.py $inputfile $version $MoE
mv $file_name /eos/user/g/gziemyte/hist_mass_point_samples
