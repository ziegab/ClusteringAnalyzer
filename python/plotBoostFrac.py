import ROOT
from ROOT import *
import glob
import numpy as np


# boost vs fraction of events passing the nClusters >= 1 cut with CLUE clusters

file_dir = '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc1p5rhoc15delt3'
print(file_dir)
root_files = glob.glob(f"{file_dir}/*.root")

# Define logarithmic bin edges
n_bins = 15  # Number of bins
x_min = 10 # Minimum x value (cannot be zero for log scale)
x_max = 1000  # Maximum x value
bin_edges = np.logspace(np.log10(x_min), np.log10(x_max), n_bins + 1)

# H_gammaPATpho = TH1F("N1", ";#gamma factor; efficiency (%)", n_bins, 0.,1000.) #numerator - gamma values of events that pass PAT pho cut in barrel
# H_gamma = TH1F("D1", ";#gamma factor; efficiency (%)", n_bins, 0.,1000.) #denominator - all gamma values of all events
H_gammaPATpho_all = TH1F("N1", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #numerator - gamma values of events that pass PAT pho cut in barrel
H_gamma_all = TH1F("D1", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #denominator - all gamma values of all events
H_gammaPATpho_all.SetStats(0)
H_gammaClust_all = TH1F("N2", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64))
H_gammaClust_all.SetStats(0)
H_gammaPATpho_EB = TH1F("N3", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #numerator - gamma values of events that pass PAT pho cut in barrel
H_gamma_EB = TH1F("D2", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #denominator - all gamma values of all events
H_gammaPATpho_EB.SetStats(0)
H_gammaClust_EB = TH1F("N4", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64))
H_gammaClust_EB.SetStats(0)
H_gammaPATpho_EE = TH1F("N5", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #numerator - gamma values of events that pass PAT pho cut in barrel
H_gamma_EE = TH1F("D3", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64)) #denominator - all gamma values of all events
H_gammaPATpho_EE.SetStats(0)
H_gammaClust_EE = TH1F("N6", ";#gamma factor; Fraction of Events", n_bins, np.array(bin_edges, dtype=np.float64))
H_gammaClust_EE.SetStats(0)

for f in root_files:#["0p1", "0p2", "0p6", "0p8", "10p0", "15p0", "20p0", "25p0", "2p0", "3p0", "4p0", "5p0", "6p0", "8p0"]:
    F = TFile(f)#TFile("/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc20delt3/hist_AtoGG_1000events"+f+"Ma2dc20rhoc3delt.root")
    T = F.Get("clus/Events")
    for e in T:
        H_gamma_all.Fill(T.gammaval[0])
        if (T.passEvent[0]>0): 
            H_gammaPATpho_all.Fill(T.gammaval[0])
        if (len(T.pEE_cluster_E)>0) or (len(T.mEE_cluster_E)>0) or (len(T.EB_cluster_E)>0):
            H_gammaClust_all.Fill(T.gammaval[0])

        if (abs(T.higgs_eta[0])<1.5):
            H_gamma_EB.Fill(T.gammaval[0])
            if (T.passEvent[0]>0):
                H_gammaPATpho_EB.Fill(T.gammaval[0])
            if (len(T.EB_cluster_E)>0):
                H_gammaClust_EB.Fill(T.gammaval[0])

        if (abs(T.higgs_eta[0])>1.5):
            H_gamma_EE.Fill(T.gammaval[0])
            if (T.passEvent[0]>0): 
                H_gammaPATpho_EE.Fill(T.gammaval[0])
            if (len(T.pEE_cluster_E)>0) or (len(T.mEE_cluster_E)>0):
                H_gammaClust_EE.Fill(T.gammaval[0])
            

H_gammaPATpho_all.Divide(H_gamma_all)
H_gammaClust_all.Divide(H_gamma_all)
H_clustoverPAT_all = H_gammaClust_all.Clone("clustoverPAT")
H_clustoverPAT_all.Divide(H_gammaPATpho_all)
H_gammaPATpho_EB.Divide(H_gamma_EB)
H_gammaClust_EB.Divide(H_gamma_EB)
H_clustoverPAT_EB = H_gammaClust_EB.Clone("clustoverPAT")
H_clustoverPAT_EB.Divide(H_gammaPATpho_EB)
H_gammaPATpho_EE.Divide(H_gamma_EE)
H_gammaClust_EE.Divide(H_gamma_EE)
H_clustoverPAT_EE = H_gammaClust_EE.Clone("clustoverPAT")
H_clustoverPAT_EE.Divide(H_gammaPATpho_EE)


C = TCanvas()
C.cd()
C.SetLogx()
H_gammaPATpho_all.SetTitle("Events with >=1 Photon (Barrel + Endcaps) - WP90 Cut")
H_gammaPATpho_all.GetYaxis().SetRangeUser(0, 1)
H_gammaPATpho_all.SetLineColor(ROOT.kBlue)
H_gammaPATpho_all.SetLineWidth(2)
H_gammaPATpho_all.Draw("hist")
C.Print("PATphoBoostdc1p5rhoc15delt3alleventsv3.root")

C10 = TCanvas()
C10.cd()
C10.SetLogx()
H_gammaPATpho_EB.SetTitle("Events with >=1 Photon (Barrel) - WP90 Cut")
H_gammaPATpho_EB.GetYaxis().SetRangeUser(0, 1)
H_gammaPATpho_EB.SetLineColor(ROOT.kBlue)
H_gammaPATpho_EB.SetLineWidth(2)
H_gammaPATpho_EB.Draw("hist")
C10.Print("PATphoBoostdc1p5rhoc15delt3EBv3.root")

C11 = TCanvas()
C11.cd()
C11.SetLogx()
H_gammaPATpho_EE.SetTitle("Events with >=1 Photon (Endcaps) - WP90 Cut")
H_gammaPATpho_EE.GetYaxis().SetRangeUser(0, 1)
H_gammaPATpho_EE.SetLineColor(ROOT.kBlue)
H_gammaPATpho_EE.SetLineWidth(2)
H_gammaPATpho_EE.Draw("hist")
C11.Print("PATphoBoostdc1p5rhoc15delt3EEv3.root")

# C2 = TCanvas()
# C2.cd()
# C2.SetLogx()
# H_gammaClust.SetTitle("Events with >=1 Photon - CLUE Clustering")
# H_gammaClust.Draw("hist")
# C2.Print("ClustBoostdc1p5rhoc15delt3v3.root")

# for H_clustoverPAT in [H_clustoverPAT_all, H_clustoverPAT_EB, H_clustoverPAT_EE]:
#     C3 = TCanvas()
#     C3.cd()
#     C3.SetLogx()
#     H_clustoverPAT.SetTitle("CLUE Clusters over PAT Photons")    
#     H_clustoverPAT.Draw("hist")
#     C3.Print("CLUEoPATdc1p5rhoc15delt3v3_{}.root")

C4 = TCanvas()
C4.cd()
C4.SetLogx()
H_gammaPATpho_all.SetLineColor(ROOT.kBlue)
H_gammaPATpho_all.SetLineWidth(2)
H_gammaClust_all.SetLineColor(ROOT.kRed)
H_gammaClust_all.SetLineWidth(2)
H_gammaPATpho_all.SetTitle("Events with >=1 Photon (Barrel + Endcaps)")
H_gammaPATpho_all.GetYaxis().SetRangeUser(0, 1)
H_gammaClust_all.GetYaxis().SetRangeUser(0, 1)
H_gammaPATpho_all.Draw("hist")
H_gammaClust_all.Draw("hist SAME")
legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)  # x1, y1, x2, y2 in normalized coordinates
legend.AddEntry(H_gammaPATpho_all, "PAT WP90 Cut", "l")  # "l" means line
legend.AddEntry(H_gammaClust_all, "CLUE Clusters", "l")
legend.Draw()
C4.Update()
C4.Print("CLUEandPATalldc1p5rhoc15delt3v1.root")

C5 = TCanvas()
C5.cd()
C5.SetLogx()
H_gammaPATpho_EB.SetLineColor(ROOT.kBlue)
H_gammaPATpho_EB.SetLineWidth(2)
H_gammaClust_EB.SetLineColor(ROOT.kRed)
H_gammaClust_EB.SetLineWidth(2)
H_gammaPATpho_EB.SetTitle("Events with >=1 Photon (Barrel)")
H_gammaPATpho_EB.GetYaxis().SetRangeUser(0, 1)
H_gammaClust_EB.GetYaxis().SetRangeUser(0, 1)
H_gammaPATpho_EB.Draw("hist")
H_gammaClust_EB.Draw("hist SAME")
legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)  # x1, y1, x2, y2 in normalized coordinates
legend.AddEntry(H_gammaPATpho_EB, "PAT WP90 Cut", "l")  # "l" means line
legend.AddEntry(H_gammaClust_EB, "CLUE Clusters", "l")
legend.Draw()
C5.Update()
C5.Print("CLUEandPATEBdc1p5rhoc15delt3v1.root")

C6 = TCanvas()
C6.cd()
C6.SetLogx()
H_gammaPATpho_EE.SetLineColor(ROOT.kBlue)
H_gammaPATpho_EE.SetLineWidth(2)
H_gammaClust_EE.SetLineColor(ROOT.kRed)
H_gammaClust_EE.SetLineWidth(2)
H_gammaPATpho_EE.SetTitle("Events with >=1 Photon (Endcaps)")
H_gammaPATpho_EE.GetYaxis().SetRangeUser(0, 1)
H_gammaClust_EE.GetYaxis().SetRangeUser(0, 1)
H_gammaPATpho_EE.Draw("hist")
H_gammaClust_EE.Draw("hist SAME")
legend = ROOT.TLegend(0.7, 0.1, 0.9, 0.3)  # x1, y1, x2, y2 in normalized coordinates
legend.AddEntry(H_gammaPATpho_EE, "PAT WP90 Cut", "l")  # "l" means line
legend.AddEntry(H_gammaClust_EE, "CLUE Clusters", "l")
legend.Draw()
C6.Update()
C6.Print("CLUEandPATEEdc1p5rhoc15delt3v1.root")