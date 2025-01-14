import ROOT
from ROOT import *
import glob
import numpy as np
import re

# delta vs. avg. number of clusters per event

file_dir = ['/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt2',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt3',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt4',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt10',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt20']

boosts = [10, 20, 50, 100, 200, 500, 700, 1000]
deltavals = [2, 3, 4, 10, 20]
histograms = []

for g in boosts:
    hname = f"boost{g}"
    htitle = f";delta value; Avg. CLUE Cluster per EB Event"
    hist = TH1F(hname, htitle, 3, 0., 21.)
    histograms.append(hist)

for g in range(len(boosts)):
    for f in range(len(file_dir)):
        root_files = glob.glob(f"{file_dir[f]}/*.root")
        EBevents = 0
        EBclusttot = 0
        for root_file in root_files:
            F = TFile(root_file)
            T = F.Get("clus/Events")
            for e in T:
                if (abs(T.higgs_eta[0])<1.5) and (boosts[g]-1) < T.gammaval[0] < (boosts[g]+1):
                    EBevents += 1
                    if (len(T.nClusters)>0):
                        EBclusttot += T.nClusters[0]
        if EBevents>0:
            if 0<deltavals[f]<5:
                histograms[g].Fill(deltavals[f], ((EBclusttot/EBevents)/3))
            else:
                histograms[g].Fill(deltavals[f], (EBclusttot/EBevents))

C = TCanvas()
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # x1, y1, x2, y2 in normalized coordinates
for idx, hist in enumerate(histograms):
    hist.SetLineColor(idx + 1)  # Different color for each
    hist.SetStats(0)
    hist.SetLineWidth(2)
    hist.SetTitle("delta values vs. CLUE Clusters")
    hist.GetYaxis().SetRangeUser(0, 2)
    legend.AddEntry(hist, f"Boost = {boosts[idx]}", "l")
    draw_option = "hist" if idx == 0 else "hist SAME"
    hist.Draw(draw_option)
legend.Draw()
C.Print("deltav1.root")
            