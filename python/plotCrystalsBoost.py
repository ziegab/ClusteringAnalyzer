import ROOT
from ROOT import *
import glob
import numpy as np


# boost vs Number of Crystals 

file_dir = '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt3'
print(file_dir)
root_files = glob.glob(f"{file_dir}/*.root")

# Define logarithmic bin edges
n_bins_x = 15  # Number of bins
n_bins_y = 9
x_min = 10 # Minimum x value (cannot be zero for log scale)
x_max = 1000  # Maximum x value
y_min = 0.1
y_max = 25
y_min_eta = -2
y_max_eta = 2
bin_edges_x = np.logspace(np.log10(x_min), np.log10(x_max), n_bins_x + 1)
bin_edges_y = np.logspace(np.log10(y_min), np.log10(y_max), n_bins_y + 1)

H_gammaMass_EB = TH2F("Boost vs. Mass Clusters", ";#gamma factor: Higgs Mass", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y, bin_edges_y)
H_gammaMasscounter_EB = TH2F("Boost vs. Mass Clusters counter", ";#gamma factor: Higgs Mass", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y, bin_edges_y)
H_gammaMass_EB.SetStats(0)
H_gammaEta_EB = TH2F("Boost vs. Eta Clusters", ";#gamma factor: Higgs Eta", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y+1, y_min_eta, y_max_eta)
H_gammaEtacounter_EB = TH2F("Boost vs. Eta Clusters counter", ";#gamma factor: Higgs Eta", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y+1, y_min_eta, y_max_eta)
H_gammaEta_EB.SetStats(0)
H_crystalplot = TH2F("Crystal Plot", "eta: phi", 170, -85, 85, 360, 0, 360)
H_crystalplot.SetStats(0)
# H_gammaMass_EE = TH2F("Boost vs. Mass Clusters", ";#gamma factor: Higgs Mass", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y, bin_edges_y)
# H_gammaMasscounter_EE = TH2F("Boost vs. Mass Clusters counter", ";#gamma factor: Higgs Mass", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y, bin_edges_y)
# H_gammaMass_EE.SetStats(0)

for f in root_files:
    F = TFile(f)
    T = F.Get("clus/Events")
    for e in T:
        if (abs(T.higgs_eta[0])<1.5):
            # Need to check if clusters are the right match (deltaR)
            if (len(T.cluster_eta)>0):
                if T.cluster_eta[0]>0:
                    deltaeta = T.higgs_eta[0] - (np.log10(np.tan(T.cluster_eta[0]/2)))
                else:
                    deltaeta = T.higgs_eta[0] + (np.log10(np.tan(abs(T.cluster_eta[0])/2)))
                deltaphi = T.higgs_phi[0] - ((T.cluster_phi[0])*(np.pi/180))
                deltaR = np.sqrt((deltaeta)**2 + (deltaphi)**2)
                if deltaR <= 0.5:
                    H_gammaMass_EB.Fill(T.gammaval[0], T.mH[0], len(T.EB_E[0]))
                    H_gammaMasscounter_EB.Fill(T.gammaval[0], T.mH[0])
                    H_gammaEta_EB.Fill(T.gammaval[0], T.higgs_eta[0], len(T.EB_E[0]))
                    H_gammaEtacounter_EB.Fill(T.gammaval[0], T.higgs_eta[0])
                    if len(T.EB_E[0]) == 133:
                        EB_eta = T.EB_eta[0]
                        EB_phi = T.EB_phi[0]
                        EB_E = T.EB_E[0]
                        for i in range(len(T.EB_E[0])):
                            H_crystalplot.Fill(EB_eta[i], EB_phi[i], EB_E[i])
                            # H_crystalplot.Fill(EB_eta[i], EB_phi[i])
        # if (abs(T.higgs_eta[0])>1.5):
        #     # Need to check if clusters are the right match (deltaR)
        #     if (len(T.cluster_eta)>0):
        #         if T.cluster_eta[0]>0:
        #             deltaeta = T.higgs_eta[0] - (np.log10(np.tan(T.cluster_eta[0]/2)))
        #         else:
        #             deltaeta = T.higgs_eta[0] + (np.log10(np.tan(abs(T.cluster_eta[0])/2)))
        #         deltaphi = T.higgs_phi[0] - ((T.cluster_phi[0])*(np.pi/180))
        #         deltaR = np.sqrt((deltaeta)**2 + (deltaphi)**2)
        #         if deltaR <= 0.5:
        #             H_gammaMass_EB.Fill(T.gammaval[0], T.mH[0], T.nClusters[0])
        #             H_gammaMasscounter_EB.Fill(T.gammaval[0], T.mH[0])

# gStyle.SetPalette(109)
# gStyle.SetPaintTextFormat("4.2f")
# C2 = TCanvas()
# C2.cd()
# C2.SetLogx()
# C2.SetLogy()
# H_gammaMasscounter_EB.SetTitle("Crystals Exist Counter (Barrel)")
# H_gammaMasscounter_EB.Draw("hist")
# H_gammaMasscounter_EB.Draw("colztext")
# xAxis = H_gammaMasscounter_EB.GetXaxis()
# xAxis.SetTitle("#gamma Factor")
# xAxis.CenterTitle(kTRUE)
# yAxis = H_gammaMasscounter_EB.GetYaxis()
# yAxis.SetTitle("Higgs Mass (GeV)")
# yAxis.CenterTitle(kTRUE)
# C2.Print("CrystalsBoostMassdc2rhoc15delt3counter.root")

H_gammaMass_EB.Divide(H_gammaMasscounter_EB)
H_gammaEta_EB.Divide(H_gammaEtacounter_EB)
gStyle.SetPalette(109)
gStyle.SetPaintTextFormat("4.2f")
C = TCanvas()
C.cd()
C.SetLogx()
C.SetLogy()
H_gammaMass_EB.SetTitle("Avg. Number of Crystals (Barrel)")
H_gammaMass_EB.Draw("hist")
H_gammaMass_EB.Draw("colztext")
xAxis = H_gammaMass_EB.GetXaxis()
xAxis.SetTitle("#gamma Factor")
xAxis.CenterTitle(kTRUE)
yAxis = H_gammaMass_EB.GetYaxis()
yAxis.SetTitle("Higgs Mass (GeV)")
yAxis.CenterTitle(kTRUE)
C.Print("CrystalsBoostMassdc2rhoc15delt3.root")

C3 = TCanvas()
C3.cd()
H_crystalplot.Draw("hist")
C3.Print("testplotcrystal1.root")


# gStyle.SetPalette(109)
# gStyle.SetPaintTextFormat("4.2f")
# C1 = TCanvas()
# C1.cd()
# C1.SetLogx()
# # C1.SetLogy()
# H_gammaEta_EB.SetTitle("Avg. Number of Crystals (Barrel)")
# H_gammaEta_EB.Draw("hist")
# H_gammaEta_EB.Draw("colztext")
# xAxis = H_gammaEta_EB.GetXaxis()
# xAxis.SetTitle("#gamma Factor")
# xAxis.CenterTitle(kTRUE)
# yAxis = H_gammaEta_EB.GetYaxis()
# yAxis.SetTitle("Higgs Eta (GeV)")
# yAxis.CenterTitle(kTRUE)
# C1.Print("CrystalsBoostEtadc2rhoc15delt3.root")