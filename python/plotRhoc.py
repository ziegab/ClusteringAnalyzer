import ROOT
from ROOT import *
import glob
import numpy as np
import re

# boost vs fraction of events passing the nClusters >= 1 cut with CLUE clusters

file_dir = ['/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc5delt3',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc10delt3',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt3',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc20delt3',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc25delt3',
            '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc30delt3']

boosts = [10, 20, 50, 100, 200, 500, 700, 1000]
rhocvals = [5, 10, 15, 20, 25, 30]
histograms = []

for g in boosts:
    hname = f"boost{g}"
    htitle = f";rho_c value; Energy Ratio"
    hist = TH1F(hname, htitle, 6, 2.5, 32.5)
    histograms.append(hist)

for g in range(len(boosts)):
    for f in range(len(file_dir)):
        root_files = glob.glob(f"{file_dir[f]}/*.root")
        totenergyratio = 0
        filecounter = 0
        # EBclusttot = 0
        for root_file in root_files:
            EBevents = 0
            # filecounter += 1
            temptotenergyratio = 0
            F = TFile(root_file)
            T = F.Get("clus/Events")
            for e in T:
                if (abs(T.higgs_eta[0])<1.5) and (boosts[g]-1) < T.gammaval[0] < (boosts[g]+1):
                    EBevents += 1
                    if (len(T.nClusters)>0):
                        EBClustEtemp = 0
                        genE = T.EH[0]
                        for i in range(T.nClusters[0]):
                            EBClustEtemp += T.cluster_E[i]
                        energyratio = EBClustEtemp/genE
                        temptotenergyratio += energyratio
            if EBevents>0:
                totenergyratio += temptotenergyratio/EBevents
                filecounter += 1
                        # EBclusttot += T.nClusters[0]
        histograms[g].Fill(rhocvals[f], (totenergyratio/filecounter))

C = TCanvas()
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)  # x1, y1, x2, y2 in normalized coordinates
for idx, hist in enumerate(histograms):
    hist.SetLineColor(idx + 1)  # Different color for each
    hist.SetStats(0)
    hist.SetLineWidth(2)
    hist.SetTitle("rho_c values vs. CLUE Cluster Energy / GEN Energy")
    hist.GetYaxis().SetRangeUser(0, 1.5)
    legend.AddEntry(hist, f"Boost = {boosts[idx]}", "l")
    draw_option = "hist" if idx == 0 else "hist SAME"
    hist.Draw(draw_option)
legend.Draw()
C.Print("rhocv1.root")
            

