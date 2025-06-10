import ROOT
from ROOT import *


# boost vs fraction of events passing the nClusters >= 1 cut with CLUE clusters


N = TH1F("N", ";#gamma factor; efficiency (%)", 50, 10.,1000.)
D = TH1F("D", ";#gamma factor; efficiency (%)", 50, 10.,1000.)
N.SetStats(0)

for f in ["0p1", "0p2", "0p6", "0p8", "10p0", "15p0", "20p0", "25p0", "2p0", "3p0", "4p0", "5p0", "6p0", "8p0"]:
    F = TFile("/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc20delt3/hist_AtoGG_1000events"+f+"Ma2dc20rhoc3delt.root")
    T = F.Get("clus/Events")
    for e in T:
        D.Fill(T.gammaval[0])
        if T.passEvent[0]>0:
            N.Fill(T.gammaval[0])

N.Divide(D)

C = TCanvas()
C.cd()
N.Draw("hist")
C.Print("passEvent.root")
    