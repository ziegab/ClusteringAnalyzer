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
            EB_eta, EB_phi, EB_E = T.EB_eta[0], T.EB_phi[0], T.EB_E[0]
            CID = T.EB_clustID[0]
            if len(T.EB_cluster_E) == 1 and T.nClusters[0] == 1 and len(genphoE) == 2:
                if genphoE[0] > 10 and genphoE[1]>10 and abs(genphoeta[0])<1.4 and abs(genphoeta[1])<1.4:
                    ERES.Fill((T.EB_cluster_E[0] - T.EH[0])/T.EH[0])
                    if (T.EB_cluster_E[0] - T.EH[0])/T.EH[0] < -0.1: continue
                    deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                    deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.EB_cluster_realeta[0],T.EB_cluster_realphi[0])
                    if deltaR0 < 0.1 and deltaR1 < 0.1:
                        n_event += 1
                        ECF1_1 = 0
                        ECF2_1 = 0
                        ECF3_1 = 0
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
                        C_2_1 = (ECF3_1 * ECF1_1) / (ECF2_1 * ECF2_1)
                        file.write(str(n_event))
                        file.write(", " + str(T.gammaval[0]))
                        file.write(", " + str(T.EH[0]))
                        file.write(", " + str(T.higgs_eta[0]))
                        file.write(", " + str(T.higgs_phi[0]))
                        file.write(", " + str(deltaR0))
                        file.write(", " + str(deltaR1))
                        file.write(", " + str(C_2_1))
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
    # C = TCanvas()
    # C.cd()
    # ERES.Draw("h")
    # C.Print("ERES_"+output+".root")
    
import sys
import glob

file_dir = str(sys.argv[1])
print(file_dir)
hist_files = glob.glob(f"{file_dir}/*.root")
file_counter = 0

for arg in hist_files:
    file_counter += 1
    PerClusterPreProc(arg, "SAIpreproc_AtoGG_5000events_" + sys.argv[2] + "_v" + str(file_counter))

# PerClusterPreProc("hist_AtoGG_1000events0p6Ma2dc15rhoc3delt.root", "test_0p6")
# PerClusterPreProc("hist_AtoGG_1000events3p0Ma2dc15rhoc3delt.root", "test_3p0")
# PerClusterPreProc("hist_AtoGG_1000events8p0Ma2dc15rhoc3delt.root", "test_8p0")