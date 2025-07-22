#!/bin/bash

histDirectory=/eos/user/g/gziemyte/hist_EE_samples/

cd /afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/hist_EE_samp/preproc
python3 /afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/python/SAI_EEpreproc.py $histDirectory 0p005_0p05