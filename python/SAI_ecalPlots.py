import ROOT
from ROOT import TFile, TH1F, TCanvas, TH2F, TMarker, gStyle, gPad
import math
import numpy as np
from itertools import combinations

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
        ieta_pho = int(ieta_pho + 0) + 0.5
    else:
        ieta_pho = int(ieta_pho + 1) + 0.5
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

def PerClusterPlot(input):
    n_passed = 0
    n_event = 0
    F = TFile(input)
    T = F.Get("clus/Events")
    for e in T:
        n_event += 1
        genphoE = T.genpho_E[0]
        genphoeta = T.genpho_eta[0]
        genphophi = T.genpho_phi[0]
        EB_clustID = T.EB_clustID[0]
        EB_eta, EB_phi, EB_E = T.EB_eta[0], T.EB_phi[0], T.EB_E[0]
        if len(T.EB_cluster_E) == 1 and T.nClusters[0] == 1 and len(genphoE) == 2:
        # if len(T.EB_cluster_E) > 1 and T.nClusters[0] > 1 and len(genphoE) == 2:
            if genphoE[0] > 10 and genphoE[1]>10 and abs(genphoeta[0])<1.4 and abs(genphoeta[1])<1.4:
                if (T.EB_cluster_E[0] - T.EH[0])/T.EH[0] < -0.1: continue
                deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])

                ECF1_1 = 0
                ECF2_1 = 0
                ECF3_1 = 0
                ECF1_2 = 0
                ECF2_2 = 0
                ECF3_2 = 0
                ECF1_3 = 0
                ECF2_3 = 0
                ECF3_3 = 0
                CID = T.EB_clustID[0]
                # for i in range(len(EB_E)-2):
                #     # ECF1 += EB_E[i]
                #     for j in range(len(EB_E)-1):
                #         # ECF2 += EB_E[i] * EB_E[j] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                #         if j == (len(EB_E)-2) and (CID[j] == CID[i] == 0):
                #             ECF1_1 += EB_E[i]
                #         for k in range(len(EB_E)):
                #                 # if CID[k] == CID[j] == CID[i] == 0:
                #                     ECF3_1 += EB_E[i] * EB_E[j] * EB_E[k] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j]) * manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k]) * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                #                     if k == (len(EB_E)-1):
                #                         if CID[i] == CID[j] == 0: ECF2_1 += EB_E[i] * EB_E[j] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                                        
                # for i in range(len(EB_E)):
                #     if CID[i] == 0:
                #         ECF1_2 += EB_E[i]
                # for i in range(len(EB_E)-1):
                #     for j in range(len(EB_E)):
                #         if CID[j] == 0 and CID[i] == 0:
                #             ECF2_2 += EB_E[i] * EB_E[j] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                # for i in range(len(EB_E)-2):
                #     for j in range(len(EB_E)-1):
                #         for k in range(len(EB_E)):
                #             if CID[k] == CID[j] == CID[i] == 0:
                #                 ECF3_2 += EB_E[i] * EB_E[j] * EB_E[k] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j]) * manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k]) * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])

                for i in range(len(EB_E)):
                    if CID[i] == 0: ECF1_3 += EB_E[i]
                    for j in range(len(EB_E)):
                        # if i != j: 
                        if CID[i] == CID[j] ==0: ECF2_3 += EB_E[i] * EB_E[j] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                        for k in range(len(EB_E)):
                            # if (i != j) and (i != k) and (j != k): 
                                if CID[i] == CID[j] == CID[k]: ECF3_3 += EB_E[i] * EB_E[j] * EB_E[k] * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j]) * manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k]) * manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                
                # ECF3_1 = 0.0
                # for i, j, k in combinations(range(len(EB_E)), 3):
                #     # pi, pj, pk = EB_E[i], EB_E[j], EB_E[k]
                    
                #     pt_prod =EB_E[i]* EB_E[j]* EB_E[k]
                    
                #     Rij = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                #     Rik = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                #     Rjk = manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k])
                    
                #     angular_prod = Rij * Rik * Rjk
                #     ECF3_1 += pt_prod * (angular_prod ** 1)

                for i in range(len(EB_E)):
                    if CID[i] == 0: 
                        ECF1_1 += EB_E[i]
                for i in range(len(EB_E)):
                    for j in range(i+1, len(EB_E)):
                        Rij = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                        if CID[i] == CID[j] ==0: 
                            ECF2_1 += EB_E[i] * EB_E[j] * Rij
                        for k in range(j+1, len(EB_E)):
                            Rik = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                            Rjk = manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k])
                            if CID[i] == CID[j] == CID[k] == 0:
                                ECF3_1 += EB_E[i] * EB_E[j] * EB_E[k] * Rij * Rik * Rjk


                print(ECF1_1, ECF2_1, ECF3_1)
                # print(ECF1_2, ECF2_2, ECF3_2)
                print(ECF1_3, ECF2_3, ECF3_3)
                C_2_1 = (ECF3_1 * ECF1_1) / (ECF2_1 * ECF2_1)
                # C_2_2 = (ECF3_2 * ECF1_2) / (ECF2_2 * ECF2_2)
                C_2_3 = (ECF3_3 * ECF1_3) / (ECF2_3 * ECF2_3)

                # if 0.5 < deltaR0 < 1 and 0.5 < deltaR1 < 1:
                # if 0.3 < deltaR0 < 0.5 and 0.3 < deltaR1 < 0.5:
                # if 0.1 < deltaR0 < 0.3 and 0.1 < deltaR1 < 0.3:
                # if deltaR0 < 0.1 and deltaR1 < 0.1:
                    # n_passed += 1
                    # if n_passed == int(sys.argv[2]):
                print(genphoeta[0],genphophi[0],genphoeta[1],genphophi[1], C_2_1, C_2_3, n_event)
                for i in range(len(T.EB_eta[0])):
                # for i in range(len(EB_clustID)):
                    # EB_eta, EB_phi, EB_E = T.EB_eta[0], T.EB_phi[0], T.EB_E[0]
                    H_ecalEB.Fill(EB_eta[i], EB_phi[i], EB_E[i])
                    # H_ecalEB.Fill(EB_eta[i], EB_phi[i], EB_clustID[i]+2)
                # # H_ecalEEp.Fill(T.pEE_x[0], T.pEE_y[0], T.pEE_E[0])
                # # H_ecalEEm.Fill(T.mEE_x[0], T.mEE_y[0], T.mEE_E[0])
                C = TCanvas()
                C.cd()
                H_ecalEB.SetTitle(f"DeltaR: {deltaR0}, {deltaR1} #gamma: {T.gammaval[0]}")
                # gStyle.SetPalette(109)
                gStyle.SetPaintTextFormat("4.2f")
                # gPad.SetLogz()
                H_ecalEB.Draw("colztext")
                pho1 = TMarker(eta_to_ieta(genphoeta[0], T.EB_cluster_eta[0], T.EB_cluster_realeta[0]), phi_to_iphi(genphophi[0], T.EB_cluster_phi[0], T.EB_cluster_realphi[0]), 29)
                pho2 = TMarker(eta_to_ieta(genphoeta[1], T.EB_cluster_eta[0], T.EB_cluster_realeta[0]), phi_to_iphi(genphophi[1], T.EB_cluster_phi[0], T.EB_cluster_realphi[0]), 29)
                clust = TMarker(T.EB_cluster_eta[0], T.EB_cluster_phi[0], 29)
                pho1.SetMarkerColor(ROOT.kRed)
                pho2.SetMarkerColor(ROOT.kRed)
                clust.SetMarkerColor(ROOT.kBlue)
                pho1.Draw()
                pho2.Draw()
                clust.Draw()
                C.Update()
                C.Print("H_ecalEB"+str(deltaR0)+"_"+str(deltaR1)+".root")
                H_ecalEB.Reset()
                # H_ecalEEm.Clear()
                # H_ecalEEp.Clear()
                break
    # print(n)
    
import sys

H_ecalEB = TH2F("H_ecalEB", "H_ecalEB", 170, -85, 85, 360, 0, 360)
H_ecalEEp = TH2F("H_ecalEEp", "H_ecalEEp", 100, 0, 100, 100, 0, 100)
H_ecalEEm = TH2F("H_ecalEEm", "H_ecalEEm", 100, 0, 100, 100, 0, 100)
H_ecalEB.SetStats(0)
H_ecalEEp.SetStats(0)
H_ecalEEm.SetStats(0)

PerClusterPlot(sys.argv[1])