import ROOT
from ROOT import TFile, TH1F, TCanvas, TH2F, TMarker, gStyle, gPad
import math
import numpy as np
from itertools import combinations

# def manualDeltaR(eta1, phi1, eta2, phi2):
#     dEta = eta1 - eta2
#     dPhi = phi1 - phi2
#     DeltaR = math.sqrt(dEta*dEta + dPhi*dPhi)
#     return DeltaR

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
    zEE = 317.0 # [m * 100]
    theta = 2 * np.arctan(np.exp(-eta))
    maxtheta = 2 * np.arctan(np.exp(-1.479))
    x = (zEE * np.tan(theta) * np.cos(phi )) * (50 / (zEE * np.tan(maxtheta))) + 49
    y = (zEE  * np.tan(theta) * np.sin(phi )) * (50 / (zEE * np.tan(maxtheta))) + 52
    return x, y



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
        # EB_clustID = T.EB_clustID[0]
        # EB_eta, EB_phi, EB_E = T.EB_eta[0], T.EB_phi[0], T.EB_E[0]
        mEE_clustID = T.mEE_clustID[0]
        mEE_x, mEE_y, mEE_E = T.mEE_x[0], T.mEE_y[0], T.mEE_E[0]
        mES1_x, mES1_y, mES1_E = T.mES1_x[0], T.mES1_y[0], T.mES1_E[0]
        mES2_x, mES2_y, mES2_E = T.mES2_x[0], T.mES2_y[0], T.mES2_E[0]
        pEE_clustID = T.pEE_clustID[0]
        pEE_x, pEE_y, pEE_E = T.pEE_x[0], T.pEE_y[0], T.pEE_E[0]
        pES1_x, pES1_y, pES1_E = T.pES1_x[0], T.pES1_y[0], T.pES1_E[0]
        pES2_x, pES2_y, pES2_E = T.pES2_x[0], T.pES2_y[0], T.pES2_E[0]
        strip = T.ESstrip[0]

        if sys.argv[4] == "mEE":
            if len(T.mEE_cluster_E) == 1 and T.nClusters[0] == 1 and len(genphoE) == 2 and len(strip) > 0:
            # if len(T.EB_cluster_E) > 1 and T.nClusters[0] > 1 and len(genphoE) == 2:
                if genphoE[0] > 10 and genphoE[1]>10 and n_event == int(sys.argv[3]):# and abs(genphoeta[0])<1.4 and abs(genphoeta[1])<1.4:
                    # if (T.EB_cluster_E[0] - T.EH[0])/T.EH[0] < -0.1: continue
                    # deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                    # deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])

                    # if 0.5 < deltaR0 < 1 and 0.5 < deltaR1 < 1:
                    # if 0.3 < deltaR0 < 0.5 and 0.3 < deltaR1 < 0.5:
                    # if 0.1 < deltaR0 < 0.3 and 0.1 < deltaR1 < 0.3:
                    # if deltaR0 < 0.1 and deltaR1 < 0.1:
                        # n_passed += 1
                        # if n_passed == int(sys.argv[2]):
                    print(genphoeta[0],genphophi[0],genphoeta[1],genphophi[1], n_event)
                    genpho0_x, genpho0_y = conversion_to_x_y(genphoeta[0], genphophi[0] + np.pi)
                    genpho1_x, genpho1_y = conversion_to_x_y(genphoeta[1], genphophi[1] + np.pi)
                    for i in range(len(mEE_clustID)):
                        H_ecalEEm.Fill(mEE_x[i], mEE_y[i], mEE_E[i])
                        H_ecalEEmclust.Fill(mEE_x[i], mEE_y[i], mEE_clustID[i]+2)
                    for i in range(len(T.mES1_x[0])):
                        # H_ecalESm1.Fill((mES1_x[i]-strip[i])/32 + 1, (mES1_y[i]-strip[i])/32 + 1, mES1_E[i])
                        H_ecalESm1.Fill((mES1_x[i]-strip[i])/32 + 1, (mES1_y[i]), mES1_E[i])
                        # print(strip[i])
                    for i in range(len(mES2_x)):
                        # H_ecalESm2.Fill((mES2_x[i]-strip[i])/32 + 1, (mES2_y[i]-strip[i])/32 + 1, mES2_E[i])
                        H_ecalESm2.Fill((mES2_x[i]), (mES2_y[i]-strip[i])/32 + 1, mES2_E[i])
                    # # H_ecalEEp.Fill(T.pEE_x[0], T.pEE_y[0], T.pEE_E[0])
                    # # H_ecalEEm.Fill(T.mEE_x[0], T.mEE_y[0], T.mEE_E[0])
                    C = TCanvas()
                    C.cd()
                    H_ecalEEm.SetTitle(f"#gamma: {T.gammaval[0]}")
                    # gStyle.SetPalette(109)
                    gStyle.SetPaintTextFormat("4.2f")
                    # gPad.SetLogz()
                    H_ecalEEm.Draw("colztext")
                    pho1 = TMarker(genpho0_x, genpho0_y, 29)
                    pho2 = TMarker(genpho1_x, genpho1_y, 29)
                    clust = TMarker(T.mEE_cluster_x[0], T.mEE_cluster_y[0], 29)
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
                    C.Print("H_ecalEEm"+sys.argv[2]+str(n_event)+".root")
                    H_ecalEEm.Reset()
                    # H_ecalEEm.Clear()
                    # H_ecalEEp.Clear()

                    C2 = TCanvas()
                    C2.cd()
                    # gStyle.SetPaintTextFormat("4.2f")
                    H_ecalESm1.Draw("colz")
                    C2.Update()
                    C2.Print("H_ecalES1m"+sys.argv[2]+str(n_event)+".root")
                    H_ecalESm1.Reset()

                    C3 = TCanvas()
                    C3.cd()
                    # gStyle.SetPaintTextFormat("4.2f")
                    H_ecalESm2.Draw("colz")
                    C3.Update()
                    C3.Print("H_ecalES2m"+sys.argv[2]+str(n_event)+".root")
                    H_ecalESm2.Reset()

                    C4 = TCanvas()
                    C4.cd()
                    gStyle.SetPaintTextFormat("4.2f")
                    H_ecalEEmclust.Draw("colztext")
                    pho1 = TMarker(genpho0_x, genpho0_y, 29)
                    pho2 = TMarker(genpho1_x, genpho1_y, 29)
                    clust = TMarker(T.mEE_cluster_x[0], T.mEE_cluster_y[0], 29)
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
                    C4.Print("H_ecalEEmclust"+sys.argv[2]+str(n_event)+".root")
                    H_ecalEEmclust.Reset()
                    break

        if sys.argv[4] == "pEE":
            if len(T.pEE_cluster_E) == 1 and T.nClusters[0] == 1 and len(genphoE) == 2 and len(strip) > 0:
            # if len(T.EB_cluster_E) > 1 and T.nClusters[0] > 1 and len(genphoE) == 2:
                if genphoE[0] > 10 and genphoE[1]>10 and n_event == int(sys.argv[3]):# and abs(genphoeta[0])<1.4 and abs(genphoeta[1])<1.4:
                    # if (T.EB_cluster_E[0] - T.EH[0])/T.EH[0] < -0.1: continue
                    # deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                    # deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])

                    # if 0.5 < deltaR0 < 1 and 0.5 < deltaR1 < 1:
                    # if 0.3 < deltaR0 < 0.5 and 0.3 < deltaR1 < 0.5:
                    # if 0.1 < deltaR0 < 0.3 and 0.1 < deltaR1 < 0.3:
                    # if deltaR0 < 0.1 and deltaR1 < 0.1:
                        # n_passed += 1
                        # if n_passed == int(sys.argv[2]):
                    print(genphoeta[0],genphophi[0],genphoeta[1],genphophi[1], n_event)
                    genpho0_x, genpho0_y = conversion_to_x_y(genphoeta[0], genphophi[0])
                    genpho1_x, genpho1_y = conversion_to_x_y(genphoeta[1], genphophi[1])
                    for i in range(len(pEE_clustID)):
                        H_ecalEEp.Fill(pEE_x[i], pEE_y[i], pEE_E[i])
                        H_ecalEEpclust.Fill(pEE_x[i], pEE_y[i], pEE_clustID[i]+2)
                    for i in range(len(T.pES1_x[0])):
                        # H_ecalESp1.Fill((pES1_x[i]-strip[i])/32 + 1, (pES1_y[i]-strip[i])/32 + 1, pES1_E[i])
                        H_ecalESp1.Fill((pES1_x[i]-strip[i])/32 + 1, (pES1_y[i]), pES1_E[i])
                        # print(strip[i])
                    for i in range(len(pES2_x)):
                        # H_ecalESp2.Fill((pES2_x[i]-strip[i])/32 + 1, (pES2_y[i]-strip[i])/32 + 1, pES2_E[i])
                        H_ecalESp2.Fill((pES2_x[i]), (pES2_y[i]-strip[i])/32 + 1, pES2_E[i])
                    # # H_ecalEEp.Fill(T.pEE_x[0], T.pEE_y[0], T.pEE_E[0])
                    # # H_ecalEEm.Fill(T.mEE_x[0], T.mEE_y[0], T.mEE_E[0])
                    C = TCanvas()
                    C.cd()
                    H_ecalEEp.SetTitle(f"#gamma: {T.gammaval[0]}")
                    # gStyle.SetPalette(109)
                    gStyle.SetPaintTextFormat("4.2f")
                    # gPad.SetLogz()
                    H_ecalEEp.Draw("colztext")
                    pho1 = TMarker(genpho0_x, genpho0_y, 29)
                    pho2 = TMarker(genpho1_x, genpho1_y, 29)
                    clust = TMarker(T.pEE_cluster_x[0], T.pEE_cluster_y[0], 29)
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
                    pho1.SetMarkerColor(ROOT.kRed)
                    pho2.SetMarkerColor(ROOT.kRed)
                    clust.SetMarkerColor(ROOT.kBlue)
                    pho1.Draw()
                    pho2.Draw()
                    clust.Draw()
                    C.Update()
                    C.Print("H_ecalEEp"+sys.argv[2]+str(n_event)+".root")
                    H_ecalEEp.Reset()
                    # H_ecalEEm.Clear()
                    # H_ecalEEp.Clear()

                    C2 = TCanvas()
                    C2.cd()
                    # gStyle.SetPaintTextFormat("4.2f")
                    H_ecalESp1.Draw("colz")
                    C2.Update()
                    C2.Print("H_ecalES1p"+sys.argv[2]+str(n_event)+".root")
                    H_ecalESp1.Reset()

                    C3 = TCanvas()
                    C3.cd()
                    # gStyle.SetPaintTextFormat("4.2f")
                    H_ecalESp2.Draw("colz")
                    C3.Update()
                    C3.Print("H_ecalES2p"+sys.argv[2]+str(n_event)+".root")
                    H_ecalESp2.Reset()

                    C4 = TCanvas()
                    C4.cd()
                    gStyle.SetPaintTextFormat("4.2f")
                    H_ecalEEpclust.Draw("colztext")
                    pho1 = TMarker(genpho0_x, genpho0_y, 29)
                    pho2 = TMarker(genpho1_x, genpho1_y, 29)
                    clust = TMarker(T.pEE_cluster_x[0], T.pEE_cluster_y[0], 29)
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
                    pho1.SetMarkerColor(ROOT.kRed)
                    pho2.SetMarkerColor(ROOT.kRed)
                    clust.SetMarkerColor(ROOT.kBlue)
                    pho1.Draw()
                    pho2.Draw()
                    clust.Draw()
                    C4.Update()
                    C4.Print("H_ecalEEpclust"+sys.argv[2]+str(n_event)+".root")
                    H_ecalEEpclust.Reset()
                    break
    # print(n)
    
import sys

H_ecalEB = TH2F("H_ecalEB", "H_ecalEB", 170, -85, 85, 360, 0, 360)
H_ecalEEm = TH2F("H_ecalEEm", "H_ecalEEm", 100, 0, 100, 100, 0, 100)
H_ecalEEmclust = TH2F("H_ecalEEmclust", "H_ecalEEmclust", 100, 0, 100, 100, 0, 100)
# H_ecalESm1 = TH2F("H_ecalESm1", "H_ecalESm1", 40, 0, 40, 1280, 0, 1280)
# H_ecalESm1 = TH2F("H_ecalESm1", "H_ecalESm1", 40, 0, 40, 40, 0, 40)
H_ecalESm1 = TH2F("H_ecalESm1", "H_ecalESm1", 40, 0, 40, 1280, 0, 1280)
# H_ecalESm2 = TH2F("H_ecalESm2", "H_ecalESm2", 40, 0, 40, 40, 0, 40)
H_ecalESm2 = TH2F("H_ecalESm2", "H_ecalESm2", 1280, 0, 1280, 40, 0, 40)
H_ecalEEp = TH2F("H_ecalEEp", "H_ecalEEp", 100, 0, 100, 100, 0, 100)
H_ecalEEpclust = TH2F("H_ecalEEpclust", "H_ecalEEpclust", 100, 0, 100, 100, 0, 100)
# H_ecalESp1 = TH2F("H_ecalESp1", "H_ecalESp1", 40, 0, 40, 1280, 0, 1280)
# H_ecalESp1 = TH2F("H_ecalESp1", "H_ecalESp1", 40, 0, 40, 40, 0, 40)
H_ecalESp1 = TH2F("H_ecalESp1", "H_ecalESp1", 40, 0, 40, 1280, 0, 1280)
# H_ecalESp2 = TH2F("H_ecalESp2", "H_ecalESp2", 40, 0, 40, 40, 0, 40)
H_ecalESp2 = TH2F("H_ecalESp2", "H_ecalESp2", 1280, 0, 1280, 40, 0, 40)
H_ecalEB.SetStats(0)
H_ecalEEp.SetStats(0)
H_ecalEEm.SetStats(0)
H_ecalESm2.SetStats(0)
H_ecalESm1.SetStats(0)
H_ecalESp2.SetStats(0)
H_ecalESp1.SetStats(0)

PerClusterPlot(sys.argv[1])