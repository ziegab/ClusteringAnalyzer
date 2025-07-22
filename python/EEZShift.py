import ROOT
from ROOT import TFile, TH1F, TCanvas, TH2F
import math

def EEZShift(eta, zPV): # ONLY APPLIES TO ENDCAPS
    z = 317. # [cm]
    theta = 2. * math.atan(math.exp(-1. * abs(eta)))
    R = z * math.tan(theta)
    if eta > 0:
        if zPV > 0:
            zp = z - abs(zPV)
        if zPV < 0:
            zp = z + abs(zPV)
        thetaprime = math.atan(R/zp)
        etaprime = -1. * math.log(math.tan(thetaprime/2.))
        return etaprime
    if eta < 0:
        if zPV > 0:
            zp = z + abs(zPV)
        if zPV < 0:
            zp = z - abs(zPV)
        thetaprime = math.atan(R/zp)
        etaprime = math.log(math.tan(thetaprime/2.))
        return etaprime


def PerClusterInfo(input):
    n_event = 0
    n_pEE = 0
    n_mEE = 0
    F = TFile(input)
    T = F.Get("clus/Events")
    for e in T:
        n_event += 1
        genphoE = T.genpho_E[0]
        genphoeta = T.genpho_eta[0]
        genphophi = T.genpho_phi[0]
        strip = T.ESstrip[0]
        PV_z = T.PV_z[0]
        if len(T.cluster_eta) > 0:
            cluster_eta = T.cluster_eta[0]
            etaprime = EEZShift(cluster_eta, PV_z)
            dEta = -1. * (cluster_eta - etaprime)
            if cluster_eta > 0:
                n_pEE += 1
                H_EEZShift_pEE.Fill(cluster_eta, PV_z, dEta)
                H_EEZShift_pEE_counts.Fill(cluster_eta, PV_z)
            if cluster_eta < 0:
                n_mEE += 1
                H_EEZShift_mEE.Fill(cluster_eta, PV_z, dEta)
                H_EEZShift_mEE_counts.Fill(cluster_eta, PV_z)
    # print("pEE events:", n_pEE)
    # print("mEE events:", n_mEE)

import sys
import glob

file_dir = str(sys.argv[1])
print(file_dir)
hist_files = glob.glob(f"{file_dir}/*.root")
file_counter = 0


H_EEZShift_pEE = TH2F("H_EEZShift_pEE", "H_EEZShift_pEE", 100, 1.5, 3.0, 100, -25., 25.)
H_EEZShift_mEE = TH2F("H_EEZShift_mEE", "H_EEZShift_mEE", 100, -3.0, -1.5, 100, -25., 25.)
H_EEZShift_pEE_counts = TH2F("H_EEZShift_pEE_counts", "H_EEZShift_pEE_counts", 100, 1.5, 3.0, 100, -25., 25.)
H_EEZShift_mEE_counts = TH2F("H_EEZShift_mEE_counts", "H_EEZShift_mEE_counts", 100, -3.0, -1.5, 100, -25., 25.)
# H_EEZShift_pEE_avg = TH2F("H_EEZShift_pEE_avg", "H_EEZShift_pEE_avg", 100, 1.5, 3.0, 100, -25., 25.)
# H_EEZShift_mEE_avg = TH2F("H_EEZShift_mEE_avg", "H_EEZShift_mEE_avg", 100, -3.0, -1.5, 100, -25., 25.)
H_EEZShift_pEE.SetStats(0)
H_EEZShift_mEE.SetStats(0)

for arg in hist_files:
    PerClusterInfo(arg)

H_EEZShift_pEE_avg = H_EEZShift_pEE.Clone("H_EEZShift_pEE")
H_EEZShift_mEE_avg = H_EEZShift_mEE.Clone("H_EEZShift_mEE")
H_EEZShift_pEE_avg.Divide(H_EEZShift_pEE_counts)
H_EEZShift_mEE_avg.Divide(H_EEZShift_mEE_counts)

# print("pEE total event number:", H_EEZShift_pEE_counts.GetEntries())
# print("mEE total event number:", H_EEZShift_mEE_counts.GetEntries())

C1 = TCanvas("C1", "C1", 800, 600)
C1.SetRightMargin(0.15)
C1.cd()
H_EEZShift_pEE_avg.GetXaxis().SetTitle("#eta of cluster w.r.t the detector origin")
H_EEZShift_pEE_avg.GetYaxis().SetTitle("z-position of the PV")
H_EEZShift_pEE_avg.GetZaxis().SetTitle("change in #eta")
H_EEZShift_pEE_avg.Draw("colz")
C1.Update()
C1.Print("H_EEZShift_pEE.root")

C2 = TCanvas("C2", "C2", 800, 600)
C2.SetRightMargin(0.15)
C2.cd()
H_EEZShift_mEE_avg.GetXaxis().SetTitle("#eta of cluster w.r.t the detector origin")
H_EEZShift_mEE_avg.GetYaxis().SetTitle("z-position of the PV")
H_EEZShift_mEE_avg.GetZaxis().SetTitle("change in #eta")
H_EEZShift_mEE_avg.Draw("colz")
C2.Update()
C2.Print("H_EEZShift_mEE.root")

C3 = TCanvas()
C3.cd()
H_EEZShift_pEE_counts.Draw("colz")
C3.Print("H_EEZShift_pEE_counts.root")

C4 = TCanvas()
C4.cd()
H_EEZShift_mEE_counts.Draw("colz")
C4.Print("H_EEZShift_mEE_counts.root")
