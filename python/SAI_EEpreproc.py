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
            strip = T.ESstrip[0]
            # pEE_eta, pEE_phi = T.pEE_eta[0], T.pEE_phi[0]
            # pES1_eta, pES1_phi = T.pES1_eta[0], T.pES1_phi[0]
            # pES2_eta, pES2_phi = T.pES2_eta[0], T.pES2_phi[0]
            # EE_x, EE_y, EE_z = T.EE_x[0], T.EE_y[0], T.EE_z[0]
            # ES_x, ES_y, ES_z = T.ES_x[0], T.ES_y[0], T.ES_z[0]
            if (len(T.mEE_cluster_E) == 1 or len(T.pEE_cluster_E) == 1) and T.nClusters[0] == 1 and len(genphoE) == 2:
                if genphoE[0] > 10 and genphoE[1]>10 and abs(genphoeta[0])>1.5 and abs(genphoeta[1])>1.5:
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
                    # ERES.Fill((T.EB_cluster_E[0] - T.EH[0])/T.EH[0])
                    # if (T.EB_cluster_E[0] - T.EH[0])/T.EH[0] < -0.1: continue
                    deltaR0 = manualDeltaR(genphoeta[0],genphophi[0],T.cluster_eta[0],T.cluster_phi[0])
                    deltaR1 = manualDeltaR(genphoeta[1],genphophi[1],T.cluster_eta[0],T.cluster_phi[0])
                    if deltaR0 < 0.1 and deltaR1 < 0.1:
                        n_event += 1
                        # ECF1_1 = 0
                        # ECF2_1 = 0
                        # ECF3_1 = 0
                        # for i in range(len(EB_E)):
                        #     if CID[i] == 0: 
                        #         ECF1_1 += EB_E[i]
                        # for i in range(len(EB_E)):
                        #     for j in range(i+1, len(EB_E)):
                        #         Rij = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[j], EB_phi[j])
                        #         if CID[i] == CID[j] ==0: 
                        #             ECF2_1 += EB_E[i] * EB_E[j] * Rij
                        #         for k in range(j+1, len(EB_E)):
                        #             Rik = manualDeltaR(EB_eta[i], EB_phi[i], EB_eta[k], EB_phi[k])
                        #             Rjk = manualDeltaR(EB_eta[j], EB_phi[j], EB_eta[k], EB_phi[k])
                        #             if CID[i] == CID[j] == CID[k] == 0:
                        #                 ECF3_1 += EB_E[i] * EB_E[j] * EB_E[k] * Rij * Rik * Rjk
                        # C_2_1 = (ECF3_1 * ECF1_1) / (ECF2_1 * ECF2_1)
                        file.write(str(n_event))
                        file.write(", " + str(T.gammaval[0]))
                        file.write(", " + str(EE_cluster_E))
                        file.write(", " + str(T.cluster_eta[0]))
                        file.write(", " + str(T.cluster_phi[0]))
                        # file.write(", " + str(deltaR0))
                        # file.write(", " + str(deltaR1))
                        # file.write(", " + str(C_2_1))
                        # file.write(", " + str(T.EB_cluster_realeta[0]))
                        # file.write(", " + str(T.EB_cluster_realphi[0]))
                        CID = EE_clustID
                        Cy = EE_y
                        Cx = EE_x
                        CE = EE_E
                        for c in range(len(CID)):
                            if CID[c] == 0:   
                                file.write(", " + str(CE[c]))
                                file.write(", " + str(int(Cx[c] - EE_cluster_x)))
                                file.write(", " + str(int(Cy[c] - EE_cluster_y)))
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
    startidx1 = arg.find("v")
    startidx2 = arg.find("v", startidx1 + 1)
    endidx = arg.rfind(".root")
    versionname = arg[(startidx2 + 1):endidx]
    PerClusterPreProc(arg, "SAIpreproc_AtoGG_5000events_" + sys.argv[2] + "_v" + str(versionname))

# PerClusterPreProc(sys.argv[1], sys.argv[2])

# PerClusterPreProc("hist_AtoGG_1000events0p6Ma2dc15rhoc3delt.root", "test_0p6")
# PerClusterPreProc("hist_AtoGG_1000events3p0Ma2dc15rhoc3delt.root", "test_3p0")
# PerClusterPreProc("hist_AtoGG_1000events8p0Ma2dc15rhoc3delt.root", "test_8p0")