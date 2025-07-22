#!/bin/bash

inputfile=
version=
MoE=
active_directory=/afs/cern.ch/user/g/gziemyte/diphogen/gitpull/ClusteringAnalyzer/
file_name="hist_AtoGG_5000events_${MoE}_Ma2dc15rhoc3deltnoresv${version}.root"

echo "CLUE Clustering"
cd $active_directory
eval `scramv1 runtime -sh`
cmsRun ../python/ConfFile_cfg.py $inputfile $version $MoE
mv $file_name /eos/user/g/gziemyte/hist_EE_samples
