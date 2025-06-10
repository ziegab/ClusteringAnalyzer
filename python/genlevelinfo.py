import ROOT
from ROOT import *
import glob
import numpy as np
import re

# taking generator level info and making plots

# file_dir = '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt3'
# root_files = glob.glob(f"{file_dir}/*.root")
root_files = glob.glob(f"/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/SAI/hist_AtoGG_1000events0p05Ma2dc15rhoc3delt.root")

H_diphoenergyfrac = TH1F("Energy Fraction", "Energy Fraction", 10, 0., 1.)
# H_diphoenergyfrac.SetStats(0)

# Define logarithmic bin edges
n_bins_x = 8  # Number of bins
n_bins_y = 5
x_min = 1 # Minimum x value (cannot be zero for log scale)
x_max = 100  # Maximum x value
y_min = 120
y_max = 130
bin_edges_x = np.logspace(np.log10(x_min), np.log10(x_max), n_bins_x + 1)
bin_edges_y = np.logspace(np.log10(y_min), np.log10(y_max), n_bins_y + 1)

n_bins_x2 = 12  # Number of bins
n_bins_y2 = 12
x2_min = 10 # Minimum x value (cannot be zero for log scale)
x2_max = 10000  # Maximum x value
y2_min = 10
y2_max = 10000
bin_edges_x2 = np.logspace(np.log10(x2_min), np.log10(x2_max), n_bins_x2 + 1)
bin_edges_y2 = np.logspace(np.log10(y2_min), np.log10(y2_max), n_bins_y2 + 1)

H_gammaMass_EB = TH2F("Boost vs. Mass Clusters", ";#gamma factor: Higgs Mass", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y, bin_edges_y)
H_gammaMasscounter_EB = TH2F("Boost vs. Mass Clusters counter", ";#gamma factor: Higgs Mass", n_bins_x, np.array(bin_edges_x, dtype=np.float64), n_bins_y, bin_edges_y)
H_gammaMass_EB.SetStats(0)

H_E1E2_EB = TH2F("E1 vs. E2", ";LargeE:SmallE", n_bins_x2, np.array(bin_edges_x2, dtype=np.float64), n_bins_y2, bin_edges_y2)
H_E1E2_EB.SetStats(0)

# ## 2D energy of each photon in diphoton pair plot
# for f in root_files:
#     F = TFile(f)
#     T = F.Get("clus/Events")
#     for e in T:
#         energyFraction = 0
#         eventPhotons = T.genpho_E[0]
#         if len(eventPhotons) == 2 and (50 <= T.gammaval[0] <= 90):
#             if eventPhotons[0] > eventPhotons[1]:
#                 H_E1E2_EB.Fill(eventPhotons[0], eventPhotons[1])
#             else:
#                 H_E1E2_EB.Fill(eventPhotons[1], eventPhotons[0])

# gStyle.SetPalette(109)
# gStyle.SetPaintTextFormat("4.2f")
# C = TCanvas()
# C.cd()
# C.SetLogx()
# C.SetLogy()
# H_E1E2_EB.SetTitle("Diphoton Energy Map (boost = 50-90)")
# H_E1E2_EB.Draw("hist")
# H_E1E2_EB.Draw("colztext")
# xAxis = H_E1E2_EB.GetXaxis()
# xAxis.SetTitle("Larger #gamma Energy (GeV)")
# xAxis.CenterTitle(kTRUE)
# yAxis = H_E1E2_EB.GetYaxis()
# yAxis.SetTitle("Smaller #gamma Energy (GeV)")
# yAxis.CenterTitle(kTRUE)
# C.Print("EnergyMapBoost90v1.root")


## 1D energy fraction plot (small energy / large energy photon)
for f in root_files:
    F = TFile(f)
    T = F.Get("clus/Events")
    for e in T:
        energyFraction = 0
        eventPhotons = T.genpho_E[0]
        if len(eventPhotons) == 2: # and (90 <= T.gammaval[0] <= 110):
            # print(eventPhotons[0], eventPhotons[1])
            if eventPhotons[0] > eventPhotons[1]:
                energyFraction = eventPhotons[1]/eventPhotons[0]
            else:
                energyFraction = eventPhotons[0]/eventPhotons[1]
            # energyFraction = abs(eventPhotons[0]-eventPhotons[1])/T.EH[0]
            H_diphoenergyfrac.Fill(energyFraction)


C = TCanvas()
C.cd()
H_diphoenergyfrac.SetTitle("Energy Fraction (mH = 0.05 GeV)")
H_diphoenergyfrac.Draw("hist")
C.Print("EnergyFractionmH0p05.root")


### 2D energy fraction plot of mass vs. boost, energy fraction is (E1-E2)/EH
# for f in root_files:
#     F = TFile(f)
#     T = F.Get("clus/Events")
#     for e in T:
#         energyFraction = 0
#         eventPhotons = T.genpho_E[0]
#         if len(eventPhotons) == 2:
#             energyFraction = abs(eventPhotons[0]-eventPhotons[1])/T.EH[0]
#             H_gammaMass_EB.Fill(T.gammaval[0], T.mH[0], energyFraction)
#             H_gammaMasscounter_EB.Fill(T.gammaval[0], T.mH[0])
        
# H_gammaMass_EB.Divide(H_gammaMasscounter_EB)
# gStyle.SetPalette(109)
# gStyle.SetPaintTextFormat("4.2f")
# C = TCanvas()
# C.cd()
# C.SetLogx()
# C.SetLogy()
# H_gammaMass_EB.SetTitle("Avg. Energy Fraction (all 2 photon events)")
# H_gammaMass_EB.Draw("hist")
# H_gammaMass_EB.Draw("colztext")
# xAxis = H_gammaMass_EB.GetXaxis()
# xAxis.SetTitle("#gamma Factor")
# xAxis.CenterTitle(kTRUE)
# yAxis = H_gammaMass_EB.GetYaxis()
# yAxis.SetTitle("Higgs Mass (GeV)")
# yAxis.CenterTitle(kTRUE)
# C.Print("BoostMassEnergyFractionAll2.root")

