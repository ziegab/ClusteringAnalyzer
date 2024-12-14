import ROOT
from ROOT import *
import glob
import numpy as np


# boost vs fraction of events passing the nClusters >= 1 cut with CLUE clusters

file_dir = '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc1p5rhoc15delt3'
print(file_dir)
root_files = glob.glob(f"{file_dir}/*.root")

# Define logarithmic bin edges
n_bins = 30  # Number of bins
x_min = 10  # Minimum x value (cannot be zero for log scale)
x_max = 1000  # Maximum x value
bin_edges = np.logspace(np.log10(x_min), np.log10(x_max), n_bins + 1)

# H_gammaPATpho = TH1F("N1", ";#gamma factor; efficiency (%)", n_bins, 0.,1000.) #numerator - gamma values of events that pass PAT pho cut in barrel
# H_gamma = TH1F("D1", ";#gamma factor; efficiency (%)", n_bins, 0.,1000.) #denominator - all gamma values of all events
H_gammaPATpho = TH1F("N1", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #numerator - gamma values of events that pass PAT pho cut in barrel
H_gamma = TH1F("D1", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #denominator - all gamma values of all events
H_gammaPATpho.SetStats(0)
H_gammaClust = TH1F("N2", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64))
H_gammaClust.SetStats(0)
# H_clustoverPAT = TH1F("N2", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64))
# H_clustoverPAT.SetStats(0)

for f in root_files:#["0p1", "0p2", "0p6", "0p8", "10p0", "15p0", "20p0", "25p0", "2p0", "3p0", "4p0", "5p0", "6p0", "8p0"]:
    F = TFile(f)#TFile("/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc20delt3/hist_AtoGG_1000events"+f+"Ma2dc20rhoc3delt.root")
    T = F.Get("clus/Events")
    for e in T:
        H_gamma.Fill(T.gammaval[0])
        if (T.passEvent[0]>0) and (-1.5<T.higgs_eta[0]<1.5):
            H_gammaPATpho.Fill(T.gammaval[0])
        if (len(T.EB_cluster_E)>0):
            H_gammaClust.Fill(T.gammaval[0])
            

H_gammaPATpho.Divide(H_gamma)
H_gammaClust.Divide(H_gamma)
H_clustoverPAT = H_gammaClust.Clone("clustoverPAT")
H_clustoverPAT.Divide(H_gammaPATpho)

C = TCanvas()
C.cd()
C.SetLogx()
H_gammaPATpho.SetTitle("Events with >=1 Photon - WP90 Cut")
H_gammaPATpho.Draw("hist")
C.Print("PATphoBoostdc1p5rhoc15delt3v1.root")

C2 = TCanvas()
C2.cd()
C2.SetLogx()
H_gammaClust.SetTitle("Events with >=1 Photon - CLUE Clustering")
H_gammaClust.Draw("hist")
C2.Print("ClustBoostdc1p5rhoc15delt3v1.root")

C3 = TCanvas()
C3.cd()
C3.SetLogx()
H_clustoverPAT.SetTitle("CLUE Clusters over PAT Photons")    
H_clustoverPAT.Draw("hist")
C3.Print("CLUEoPATdc1p5rhoc15delt3v1.root")