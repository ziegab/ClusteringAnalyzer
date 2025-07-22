#!/bin/bash

inputfile=/eos/user/g/gziemyte/DiphotonSamplesNoResonance/EEsamples/AtoGG_5000events0p005_0p02MoE10p0_200p0EPileupNoResv1_MiniAODpEE.root
version=pEE1
MoE=0p005_0p02
active_directory=/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/hist_EE_samp
file_name="hist_AtoGG_5000events_${MoE}_Ma2dc15rhoc3deltnoresv${version}.root"

echo "CLUE Clustering"
cd $active_directory
eval `scramv1 runtime -sh`
cmsRun ../python/ConfFile_cfg.py $inputfile $version $MoE
mv $file_name /eos/user/g/gziemyte/hist_EE_samples
