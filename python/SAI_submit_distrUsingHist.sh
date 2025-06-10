#!/bin/bash

histDirectory=/eos/user/g/gziemyte/histdc2rhoc15delt3/MoE_0p01_0p05_v2/test/

cd /afs/cern.ch/user/g/gziemyte/public
python3 SAI_distrUsingHist.py $histDirectory
