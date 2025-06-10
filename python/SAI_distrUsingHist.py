import ROOT
from ROOT import TFile, TH1F, TCanvas, TH2F, TMarker, gStyle, gPad
import math
import numpy as np
import glob

def manualDeltaR(eta1, phi1, eta2, phi2):
    dEta = eta1 - eta2
    dPhi = phi1 - phi2
    DeltaR = math.sqrt(dEta*dEta + dPhi*dPhi)
    return DeltaR

def eta_to_ieta(eta_pho, ieta_clust, eta_clust):
    ieta_pho = ieta_clust * (eta_pho/eta_clust)
    # if 0 < ieta_pho < 28:
    #     ieta_pho = int(ieta_pho + 1) + 0.5
    # elif 28 < ieta_pho < 86:
    #     ieta_pho = int(ieta_pho - 1) + 0.5
    if ieta_pho > 0:
        ieta_pho = int(ieta_pho + 3) + 0.5
    else:
        ieta_pho = int(ieta_pho - 3) + 0.5
    return ieta_pho
    # ieta = abs(eta_pho)/0.0175 + 0.5
    # if eta_pho < 0:
    #     ieta = -ieta
    # return int(ieta)

def phi_to_iphi(phi_pho, iphi_clust, phi_clust):
    # iphi_pho = iphi_clust * (phi_pho/phi_clust)
    # return (iphi_pho + 0.5)
    phi_pho_rad = phi_pho + np.pi
    phi_pho_deg = (phi_pho/np.pi)
    iphi = phi_pho_deg * 180
    if phi_pho < 0:
        iphi = iphi + 360
    iphi = iphi + 11
    return (int(iphi) + 0.5)

    
import sys

H_deltaR0 = TH1F("H_deltaR0", "DeltaR for Photon 1", 100, 0.0, 0.5)
H_deltaR1 = TH1F("H_deltaR1", "DeltaR for Photon 2", 100, 0.0, 0.5)
H_boost = TH1F("H_boost", "Boost Distribution", 50, 0, 120)
H_C2 = TH1F("H_C2", "$C_2$ Distribution", 50, 0, 10)

file_dir = str(sys.argv[1])
print(file_dir)
csv_files = glob.glob(f"{file_dir}/*.root")
print(len(csv_files))

for arg in csv_files:
    n_passed = 0
    n_event = 0
    F = TFile(arg)
    T = F.Get("clus/Events")
    for e in T:
        n_event += 1
        genphoE = T.genpho_E[0]
        genphoeta = T.genpho_eta[0]
        genphophi = T.genpho_phi[0]
        EB_clustID = T.EB_clustID[0]
        EB_eta, EB_phi, EB_E = T.EB_eta[0], T.EB_phi[0], T.EB_E[0]
        if len(T.EB_cluster_E) == 1 and T.nClusters[0] == 1 and len(genphoE) == 2:
            if genphoE[0] > 10 and genphoE[1]>10 and abs(genphoeta[0])<1.4 and abs(genphoeta[1])<1.4:
                if (T.EB_cluster_E[0] - T.EH[0])/T.EH[0] < -0.1: continue
                deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                gam = T.gammaval[0]
                n_passed += 1

                ECF1 = 0
                ECF2 = 0
                ECF3 = 0
                CID = T.EB_clustID[0]
                # for i in range(len(EB_E)):
                #     if CID[i] == 0: ECF1 += EB_E[i]
                #     for j in range(len(EB_E)):
                #         # if i != j: 
                #         if CID[i] == CID[j] == 0: ECF2 += EB_E[i] * EB_E[j] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                #         for k in range(len(EB_E)):
                #             # if (i != j) and (i != k) and (j != k): 
                #                 if CID[i] == CID[j] == CID[k] == 0: ECF3 += EB_E[i] * EB_E[j] * EB_E[k] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j]) * manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k]) * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                # for i in range(len(EB_E)):
                #     if CID[i] == 0:
                #         ECF1 += EB_E[i]
                # for i in range(len(EB_E)-1):
                #     for j in range(len(EB_E)):
                #         if CID[j] == 0 and CID[i] == 0:
                #             ECF2 += EB_E[i] * EB_E[j] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                # for i in range(len(EB_E)-2):
                #     for j in range(len(EB_E)-1):
                #         for k in range(len(EB_E)):
                #             if CID[k] == CID[j] == CID[i] == 0:
                #                 ECF3 += EB_E[i] * EB_E[j] * EB_E[k] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j]) * manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k]) * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                for i in range(len(EB_E)):
                    if CID[i] == 0: 
                        ECF1 += EB_E[i]
                for i in range(len(EB_E)):
                    for j in range(i+1, len(EB_E)):
                        Rij = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                        if CID[i] == CID[j] ==0: 
                            ECF2 += EB_E[i] * EB_E[j] * Rij
                        for k in range(j+1, len(EB_E)):
                            Rik = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                            Rjk = manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k])
                            if CID[i] == CID[j] == CID[k] == 0:
                                ECF3 += EB_E[i] * EB_E[j] * EB_E[k] * Rij * Rik * Rjk

                # print(ECF1, ECF2, ECF3)
                C_2 = (ECF3 * ECF1) / (ECF2 * ECF2)

                H_deltaR0.Fill(deltaR0)
                H_deltaR1.Fill(deltaR0)
                H_boost.Fill(gam)
                H_C2.Fill(C_2)
                # if n_passed == 500:
                #     break
                

C1 = TCanvas()
C1.cd()
H_deltaR0.Draw("h")
C1.Print("H_deltaR0.root")

C2 = TCanvas()
C2.cd()
H_deltaR1.Draw("h")
C2.Print("H_deltaR1.root")

C3 = TCanvas()
C3.cd()
H_boost.Draw("h")
C3.Print("H_boost.root")

C4 = TCanvas()
C4.cd()
H_C2.Draw("h")
C4.Print("H_C2.root")