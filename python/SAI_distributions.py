import ROOT
from ROOT import TH2F, TCanvas, TH1F
import numpy as np
import glob

def getinfo(l):
    gam = (float(l[1]))
    E = (float(l[2]))
    eta = float(l[3])
    phi = float(l[4])
    deltaR0 = float(l[5])
    deltaR1 = float(l[6])
    n_bins = 32
    length = float(int(n_bins/2))
    mH = E/gam
    # H = TH2F("H", "#gamma = "+str(gam)+", E = "+str(E)+ " GeV, #eta = "+"{:.2f}".format(round(eta, 2))+", #phi = "+"{:.2f}".format(round(phi, 2))+";#eta;#phi", n_bins,-length,length,n_bins,-length,length)
    # H.SetStats(0)
    # for i in range(6, len(l), 3):
    #     H.Fill(int(l[i+1]),int(l[i+2]),float(l[i]))
    # H_deltaR0 = TH1F("H_deltaR0", "DeltaR for Photon 1", 100, 0.0, 0.5)
    # H_deltaR1 = TH1F("H_deltaR1", "DeltaR for Photon 2", 100, 0.0, 0.5)
    # H_boost = TH1F("H_boost", "Boost Distribution", 100, 0, 230)
    # H_MoE = TH1F("H_MoE", "MoE Distribution", 100, 0, 1.5)

    # H_deltaR0.Fill(deltaR0)
    # H_deltaR1.Fill(deltaR1)
    # H_boost.Fill(gam)
    # H_MoE.Fill(1/gam)
    return gam, deltaR0, deltaR1

import csv
import sys

H_deltaR0 = TH1F("H_deltaR0", "DeltaR for Photon 1", 100, 0.0, 0.5)
H_deltaR1 = TH1F("H_deltaR1", "DeltaR for Photon 2", 100, 0.0, 0.5)
H_boost = TH1F("H_boost", "Boost Distribution", 50, 0, 120)
H_MoE = TH1F("H_MoE", "MoE Distribution", 50, 0, 0.1)

file_dir = str(sys.argv[1])
print(file_dir)
csv_files = glob.glob(f"{file_dir}/*.csv")
print(len(csv_files))

for arg in csv_files:
    with open(arg) as csvfile:
        reader = csv.reader(csvfile)
        # mH = float(get_mass(str(sys.argv[1]), "SAIpreproc_AtoGG_"))
        # output = sys.argv[2]
        n = 0
        for row in reader:
            n+=1
            # makeentry(row, output, mH)
            gam, deltaR0, deltaR1 = getinfo(row)

            if deltaR0 < 0.1 and deltaR1 < 0.1:
                H_deltaR0.Fill(deltaR0)
                H_deltaR1.Fill(deltaR1)
                H_boost.Fill(gam)
                H_MoE.Fill(1/gam)

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
H_MoE.Draw("h")
C4.Print("H_MoE.root")