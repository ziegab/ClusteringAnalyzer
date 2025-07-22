#!/bin/bash

inputfile=/eos/user/g/gziemyte/AtoGG_5000events0p1Ma14p0_18p0EPileupNoResv1_MiniAODvalidationSample.root
version=14p0_18p0E
MoE=0p1Ma
active_directory=/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/hist_mass_point
file_name="hist_AtoGG_5000events_${MoE}_Ma2dc15rhoc3deltnoresv${version}.root"

echo "CLUE Clustering"
cd $active_directory
eval `scramv1 runtime -sh`
cmsRun ../python/ConfFile_cfg.py $inputfile $version $MoE
mv $file_name /eos/user/g/gziemyte/hist_mass_point_samples
