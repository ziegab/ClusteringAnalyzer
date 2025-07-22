import ROOT
from ROOT import TH2F, TCanvas
import numpy as np

def makeplot(l, n):
    index = int(l[0])
    gam = int(float(l[1]))
    E = int(float(l[2]))
    eta = float(l[3])
    phi = float(l[4])
    # mass = float(l[1])/float(l[0])
    n_bins = 32
    length = float(int(n_bins/2))
    H = TH2F("H", "#gamma = "+str(gam)+", E = "+str(E)+ " GeV, #eta = "+"{:.2f}".format(round(eta, 2))+", #phi = "+"{:.2f}".format(round(phi, 2))+";#eta;#phi", n_bins,-length,length,n_bins,-length,length)
    H.SetStats(0)
    for i in range(5, len(l), 3):
        H.Fill(int(l[i+1]),int(l[i+2]),float(l[i]))
    # maxcrystalE = H.GetBinContent(8, 8)
    if float(l[1]) > 0:
        H.Scale(1.0 / float(l[2]))

    C = TCanvas()
    C.cd()
    C.SetLogz()
    H.Draw("col")
    csvfilename = str(sys.argv[1])
    csvfilename = csvfilename.rstrip(".csv")
    C.Print(csvfilename+str(n)+".root")
    # array = H_rotated.GetArray()
    # for y in range(n_bins,0,-1):
    #     # for x in range(17,0,-1):
    #     for x in range(1, (n_bins+1)):
    #         index = x  + y * (n_bins+2)
    #         print(f"{array[index]:.2f}", end=" ")
    #     print("")

import csv
import sys

with open(sys.argv[1]) as csvfile:
    reader = csv.reader(csvfile)
    n = 0
    for row in reader:
        n+=1
        # if 31 < int(float(row[0])) < 41:
        #     print(n)
        if n == int(sys.argv[2]):
            makeplot(row, n)