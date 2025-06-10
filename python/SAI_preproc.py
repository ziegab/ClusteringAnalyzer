import ROOT
from ROOT import TFile, TH1F, TCanvas
import math

def manualDeltaR(eta1, phi1, eta2, phi2):
    dEta = eta1 - eta2
    dPhi = phi1 - phi2
    DeltaR = math.sqrt(dEta*dEta + dPhi*dPhi)
    return DeltaR

def PerClusterPreProc(input, output):
    n_event = 0
    ERES = TH1F("ERES_"+output, ";#frac{E_{rec} - E_{gen}}{E_{gen}}", 50, -1.,1.)
    F = TFile(input)
    T = F.Get("clus/Events")
    with open(output+'.csv', mode='w') as file:
        for e in T:
            genphoE = T.genpho_E[0]
            genphoeta = T.genpho_eta[0]
            genphophi = T.genpho_phi[0]
            if len(T.EB_cluster_E) == 1 and T.nClusters[0] == 1 and len(genphoE) == 2:
                if genphoE[0] > 10 and genphoE[1]>10 and abs(genphoeta[0])<1.4 and abs(genphoeta[1])<1.4:
                    ERES.Fill((T.EB_cluster_E[0] - T.EH[0])/T.EH[0])
                    if (T.EB_cluster_E[0] - T.EH[0])/T.EH[0] < -0.1: continue
                    deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                    deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                    if deltaR0 < 0.1 and deltaR1 < 0.1:
                        n_event += 1
                        file.write(str(n_event))
                        file.write(", " + str(T.gammaval[0]))
                        file.write(", " + str(T.EH[0]))
                        file.write(", " + str(T.higgs_eta[0]))
                        file.write(", " + str(T.higgs_phi[0]))
                        file.write(", " + str(deltaR0))
                        file.write(", " + str(deltaR1))
                        # file.write(", " + str(T.EB_cluster_realeta[0]))
                        # file.write(", " + str(T.EB_cluster_realphi[0]))
                        CID = T.EB_clustID[0]
                        Cphi = T.EB_phi[0]
                        Ceta = T.EB_eta[0]
                        CE = T.EB_E[0]
                        for c in range(len(CID)):
                            if CID[c] == 0:   
                                file.write(", " + str(CE[c]))
                                file.write(", " + str(int(Ceta[c] - T.EB_cluster_eta[0])))
                                file.write(", " + str(int(Cphi[c] - T.EB_cluster_phi[0])))
                        file.write("\n")
    C = TCanvas()
    C.cd()
    ERES.Draw("h")
    C.Print("ERES_"+output+".root")
    
import sys

PerClusterPreProc(sys.argv[1], sys.argv[2])

# PerClusterPreProc("hist_AtoGG_1000events0p6Ma2dc15rhoc3delt.root", "test_0p6")
# PerClusterPreProc("hist_AtoGG_1000events3p0Ma2dc15rhoc3delt.root", "test_3p0")
# PerClusterPreProc("hist_AtoGG_1000events8p0Ma2dc15rhoc3delt.root", "test_8p0")