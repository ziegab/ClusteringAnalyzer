import ROOT
from ROOT import TH2F, TCanvas
import numpy as np

# take histograms of each cluster image (found in preproc csv file) and just put on one histogram for a given mass and boost range!
# boost = 0-300, 300-600, 600-900, 900-1200
# this script should produce 4 histogram images for a given mass point

def add_to_stacked_hist(l, n, n1, n2, n3, n4):
    gam = int(float(l[0]))
    E = float(l[1])
    eta = float(l[2])
    phi = float(l[3])

    for i in range(4, len(l), 3):
        if 0 <= gam < 250:
            n1 += 1
            if n1 < 500:
                hists_0_300[n].Fill(int(l[i+1]),int(l[i+2]),(float(l[i])/E))
        if 250 <= gam < 500:
            n2 += 1
            if n2 < 500:
                hists_300_600[n].Fill(int(l[i+1]),int(l[i+2]),(float(l[i])/E))
        if 500 <= gam < 750:
            n3 += 1
            if n3 < 500:
                hists_600_900[n].Fill(int(l[i+1]),int(l[i+2]),(float(l[i])/E))
        if 750 <= gam < 1000:
            n4 += 1
            if n4 < 500:
                hists_900_1200[n].Fill(int(l[i+1]),int(l[i+2]),(float(l[i])/E))

import csv
import sys

csvfiles = [sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]]
mH = [1.0, 2.0, 3.0, 4.0]

n_bins = 32
length = float(int(n_bins/2))
hists_0_300 = [TH2F(f"H_m{i}_0_300", f"#gamma = 0-250, mH = {str(mH[i])}"+";#eta;#phi", n_bins,-length,length,n_bins,-length,length) for i in range(4)]
hists_300_600 = [TH2F(f"H_m{i}_300_600", f"#gamma = 250-500, mH = {str(mH[i])}"+";#eta;#phi", n_bins,-length,length,n_bins,-length,length) for i in range(4)]
hists_600_900 = [TH2F(f"H_m{i}_600_900", f"#gamma = 500-750, mH = {str(mH[i])}"+";#eta;#phi", n_bins,-length,length,n_bins,-length,length) for i in range(4)]
hists_900_1200 = [TH2F(f"H_m{i}_900_1200", f"#gamma = 750-1000, mH = {str(mH[i])}"+";#eta;#phi", n_bins,-length,length,n_bins,-length,length) for i in range(4)]

for i, csvitem in enumerate(csvfiles):
    n1 = 0
    n2 = 0
    n3 = 0
    n4 = 0
    with open(csvitem) as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            add_to_stacked_hist(row, i, n1, n2, n3, n4)


# csvfilename = str(sys.argv[1])
# csvfilename = csvfilename.rstrip(".csv")

# canvas = ROOT.TCanvas("canvas", csvfilename+" Stacked Boost Plots", 1200, 600)
# canvas.Divide(4, 1)  # 4 columns, 1 row
# canvas.cd(1)
# ROOT.gPad.SetLogz(1)
# H_0_300.Draw("col")
# canvas.cd(2)
# ROOT.gPad.SetLogz(2)
# H_300_600.Draw("col")
# canvas.cd(3)
# ROOT.gPad.SetLogz(3)
# H_600_900.Draw("col")
# canvas.cd(4)
# ROOT.gPad.SetLogz(4)
# H_900_1200.Draw("col")

# canvas.Print(csvfilename+"stackboost.root")

canvas = ROOT.TCanvas("canvas", "Stacked Boost Plots", 1200, 1200)
canvas.Divide(4, 4)  # 4 columns, 4 row

for i in range(4):
        canvas.cd((i*4)+1), ROOT.gPad.SetLogz((i*4)+1)
        hists_0_300[i].SetStats(0)
        hists_0_300[i].SetBinContent(17, 17, 0)
        hists_0_300[i].Draw("col")
        canvas.cd((i*4)+2), ROOT.gPad.SetLogz((i*4)+2)
        hists_300_600[i].SetStats(0)
        hists_300_600[i].SetBinContent(17, 17, 0)
        hists_300_600[i].Draw("col")
        canvas.cd((i*4)+3), ROOT.gPad.SetLogz((i*4)+3)
        hists_600_900[i].SetStats(0)
        hists_600_900[i].SetBinContent(17, 17, 0)
        hists_600_900[i].Draw("col")
        canvas.cd((i*4)+4), ROOT.gPad.SetLogz((i*4)+4)
        hists_900_1200[i].SetStats(0)
        hists_900_1200[i].SetBinContent(17, 17, 0)
        hists_900_1200[i].Draw("col")

canvas.Print("StackBoosts.root")
