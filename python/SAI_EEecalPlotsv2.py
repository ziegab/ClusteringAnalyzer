import ROOT
from ROOT import TFile, TH1F, TCanvas, TH2F, TMarker, gStyle, gPad
import math
import numpy as np
from itertools import combinations
from scipy.ndimage import zoom

def manualDeltaR(eta1, phi1, eta2, phi2):
    dEta = eta1 - eta2
    dPhi = phi1 - phi2
    DeltaR = math.sqrt(dEta*dEta + dPhi*dPhi)
    return DeltaR

# def eta_to_ieta(eta_pho, ieta_clust, eta_clust):
#     ieta_pho = ieta_clust * (eta_pho/eta_clust)
#     # if 0 < ieta_pho < 28:
#     #     ieta_pho = int(ieta_pho + 1) + 0.5
#     # elif 28 < ieta_pho < 86:
#     #     ieta_pho = int(ieta_pho - 1) + 0.5
#     if ieta_pho > 0:
#         ieta_pho = int(ieta_pho + 0) + 0.5
#     else:
#         ieta_pho = int(ieta_pho + 1) + 0.5
#     return ieta_pho
#     # ieta = abs(eta_pho)/0.0175 + 0.5
#     # if eta_pho < 0:
#     #     ieta = -ieta
#     # return int(ieta)

def conversion_to_x_y(eta, phi):
    zEE = 317.0 # [cm]
    theta = 2 * np.arctan(np.exp(-eta))
    maxtheta = 2 * np.arctan(np.exp(-1.479))
    x = (zEE * np.tan(theta) * np.cos(phi )) * (50 / (zEE * np.tan(maxtheta))) + 50
    y = (zEE  * np.tan(theta) * np.sin(phi )) * (50 / (zEE * np.tan(maxtheta))) + 50
    return x, y

def projection_to_ES(x, y, z, ESz):
    # v = np.array([x, y, z])
    # v_unit = v / np.linalg.norm(v)
    # # ESz = 303.5
    # scale = ESz/v_unit[2]
    # xES = v_unit[0] * scale
    # yES = v_unit[1] * scale
    theta = np.arctan(y/z)
    yES = ESz * np.tan(theta)
    phi = np.arctan(x/z)
    xES = ESz * np.tan(phi)
    return xES, yES


# def phi_to_iphi(phi_pho, iphi_clust, phi_clust):
#     # iphi_pho = iphi_clust * (phi_pho/phi_clust)
#     # return (iphi_pho + 0.5)
#     phi_pho_rad = phi_pho + np.pi
#     phi_pho_deg = (phi_pho/np.pi)
#     iphi = phi_pho_deg * 180
#     if phi_pho < 0:
#         iphi = iphi + 360
#     iphi = iphi + 11
#     return (int(iphi) + 0.5)

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
        strip = T.ESstrip[0]

        if (len(T.mEE_cluster_E) == 1 or len(T.pEE_cluster_E) == 1) and T.nClusters[0] == 1 and len(genphoE) == 2: #and len(strip) > 0:
        # if len(T.EB_cluster_E) > 1 and T.nClusters[0] > 1 and len(genphoE) == 2:
            if genphoE[0] > 10 and genphoE[1]>10 and abs(genphoeta[0])>1.4 and abs(genphoeta[1])>1.4:#and n_event == int(sys.argv[3]):# and abs(genphoeta[0])<1.4 and abs(genphoeta[1])<1.4:
                # n_passed += 1
                if genphoeta[0] < 0 and len(T.pEE_cluster_E) == 0:
                    EE_clustID = T.mEE_clustID[0]
                    EE_x, EE_y, EE_E = T.mEE_x[0], T.mEE_y[0], T.mEE_E[0]
                    ES1_x, ES1_y, ES1_E = T.mES1_x[0], T.mES1_y[0], T.mES1_E[0]
                    ES2_x, ES2_y, ES2_E = T.mES2_x[0], T.mES2_y[0], T.mES2_E[0]
                elif genphoeta[0] > 0 and len(T.mEE_cluster_E) == 0:
                    EE_clustID = T.pEE_clustID[0]
                    EE_x, EE_y, EE_E = T.pEE_x[0], T.pEE_y[0], T.pEE_E[0]
                    ES1_x, ES1_y, ES1_E = T.pES1_x[0], T.pES1_y[0], T.pES1_E[0]
                    ES2_x, ES2_y, ES2_E = T.pES2_x[0], T.pES2_y[0], T.pES2_E[0]
                else:
                    continue
                EE_cluster_x, EE_cluster_y, EE_cluster_E = T.cluster_x[0], T.cluster_y[0], T.cluster_E[0]
                genpho0_x, genpho0_y = conversion_to_x_y(genphoeta[0], genphophi[0] + np.pi)
                genpho1_x, genpho1_y = conversion_to_x_y(genphoeta[1], genphophi[1] + np.pi)
                deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.cluster_eta[0],T.cluster_phi[0])
                deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.cluster_eta[0],T.cluster_phi[0])
                deltaR0_xy = manualDeltaR(genpho0_x,genpho0_y,EE_cluster_x,EE_cluster_y)
                deltaR1_xy = manualDeltaR(genpho1_x,genpho1_y,EE_cluster_x,EE_cluster_y)
                if deltaR0 > 0.5 and deltaR1 > 0.5: # and T.gammaval[0] > 60:
                # if 0.1 < deltaR0 < 0.5 and 0.1 < deltaR1 < 0.5: # and T.gammaval[0] > 60:
                    n_passed += 1 
                    if n_passed == int(sys.argv[3]):
                        print(genphoeta[0],genphophi[0],genphoeta[1],genphophi[1], n_event)
                        print(T.cluster_eta[0], T.cluster_phi[0])
                        print(deltaR0, deltaR1, deltaR0_xy, deltaR1_xy)
                        for i in range(len(EE_clustID)):
                            H_ecalEE.Fill(EE_x[i], EE_y[i], EE_E[i])
                            H_ecalEEclust.Fill(EE_x[i], EE_y[i], EE_clustID[i]+2)
                        for i in range(len(ES1_x)):
                            H_ecalES1.Fill((ES1_x[i]-strip[i])/32 + 1, (ES1_y[i]), ES1_E[i])
                        for i in range(len(ES2_x)):
                            H_ecalES2.Fill((ES2_x[i]), (ES2_y[i]-strip[i])/32 + 1, ES2_E[i])
                        
                        C = TCanvas()
                        C.cd()
                        H_ecalEE.SetTitle(f"#gamma: {T.gammaval[0]}")
                        # gStyle.SetPalette(109)
                        gStyle.SetPaintTextFormat("4.2f")
                        # gPad.SetLogz()
                        H_ecalEE.Draw("colztext")
                        pho1 = TMarker(genpho0_x, genpho0_y, 29)
                        pho2 = TMarker(genpho1_x, genpho1_y, 29)
                        clust = TMarker(EE_cluster_x, EE_cluster_y, 29)
                        pho1.SetMarkerColor(ROOT.kRed)
                        pho2.SetMarkerColor(ROOT.kRed)
                        clust.SetMarkerColor(ROOT.kBlue)
                        pho1.Draw()
                        pho2.Draw()
                        clust.Draw()
                        innerradlength = 317*np.tan(2 * np.arctan(np.exp(-3.0)))*(50 / (317 * np.tan(2 * np.arctan(np.exp(-1.479)))))
                        innerrad = ROOT.TEllipse(50, 50, innerradlength, innerradlength)  
                        innerrad.SetLineColor(ROOT.kBlack)
                        innerrad.SetLineWidth(2)
                        innerrad.SetFillStyle(0)  # transparent
                        innerrad.Draw("same")
                        outerradlength = 50
                        outerrad = ROOT.TEllipse(50, 50, outerradlength, outerradlength)  
                        outerrad.SetLineColor(ROOT.kBlack)
                        outerrad.SetLineWidth(2)
                        outerrad.SetFillStyle(0)  # transparent
                        outerrad.Draw("same")
                        C.Update()
                        C.Print("H_ecalEE"+sys.argv[2]+str(n_event)+".root")
                        H_ecalEE.Reset()

                        # C2 = TCanvas()
                        # C2.cd()
                        # # gStyle.SetPaintTextFormat("4.2f")
                        # H_ecalES1.Draw("colz")
                        # C2.Update()
                        # C2.Print("H_ecalES1m"+sys.argv[2]+str(n_event)+".root")
                        # H_ecalES1.Reset()

                        # C3 = TCanvas()
                        # C3.cd()
                        # # gStyle.SetPaintTextFormat("4.2f")
                        # H_ecalES2.Draw("colz")
                        # C3.Update()
                        # C3.Print("H_ecalES2m"+sys.argv[2]+str(n_event)+".root")
                        # H_ecalES2.Reset()

                        C4 = TCanvas()
                        C4.cd()
                        gStyle.SetPaintTextFormat("4.2f")
                        H_ecalEEclust.Draw("colztext")
                        pho1 = TMarker(genpho0_x, genpho0_y, 29)
                        pho2 = TMarker(genpho1_x, genpho1_y, 29)
                        clust = TMarker(EE_cluster_x, EE_cluster_y, 29)
                        pho1.SetMarkerColor(ROOT.kRed)
                        pho2.SetMarkerColor(ROOT.kRed)
                        clust.SetMarkerColor(ROOT.kBlue)
                        pho1.Draw()
                        pho2.Draw()
                        clust.Draw()
                        innerradlength = 317*np.tan(2 * np.arctan(np.exp(-3.0)))*(50 / (317 * np.tan(2 * np.arctan(np.exp(-1.479)))))
                        innerrad = ROOT.TEllipse(50, 50, innerradlength, innerradlength)  
                        innerrad.SetLineColor(ROOT.kBlack)
                        innerrad.SetLineWidth(2)
                        innerrad.SetFillStyle(0)  # transparent
                        innerrad.Draw("same")
                        outerradlength = 50
                        outerrad = ROOT.TEllipse(50, 50, outerradlength, outerradlength)  
                        outerrad.SetLineColor(ROOT.kBlack)
                        outerrad.SetLineWidth(2)
                        outerrad.SetFillStyle(0)  # transparent
                        outerrad.Draw("same")
                        C4.Update()
                        C4.Print("H_ecalEEclust"+sys.argv[2]+str(n_event)+".root")
                        H_ecalEEclust.Reset()
                        break

import sys

H_ecalEE = TH2F("H_ecalEE", "H_ecalEE", 100, 0, 100, 100, 0, 100)
H_ecalEEclust = TH2F("H_ecalEEclust", "H_ecalEEclust", 100, 0, 100, 100, 0, 100)
H_ecalES1 = TH2F("H_ecalES1", "H_ecalES1", 40, 0, 40, 1280, 0, 1280)
H_ecalES2 = TH2F("H_ecalES2", "H_ecalES2", 1280, 0, 1280, 40, 0, 40)
H_ecalEE.SetStats(0)
H_ecalEEclust.SetStats(0)
H_ecalES2.SetStats(0)
H_ecalES1.SetStats(0)

PerClusterPlot(sys.argv[1])