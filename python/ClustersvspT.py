# Standard Includes: ----------------------------------------#
import os
import array
from array import *
import glob
import math
from math import *
import ROOT
from ROOT import *
import sys
import csv
import itertools
from itertools import *
from optparse import OptionParser
# Obviously can only be run in a CMSSW Framework
from DataFormats.FWLite import * 
from HLTrigger import *
from ctypes import POINTER, c_int
import matplotlib.pyplot as plt
import numpy as np
# -----------------------------------------------------------#

# Goal: Take cluster and pt info from csv and plot it

# Command line: python3 python/ClustersvspT.py clusters.csv clusters.pdf 0 1000
# arguments: Name of csv file, Name of pdf file, Emin, Emax

# Initialize 1D histograms
# H_clusters = TH1F("ClustersvspT", "Higgs pT;Number of PAT photons", 100, 0, 1000)

# Pull from csv and put into lists
rows = []
numpatpho = []
patphopt = []
patphoM = []
patphoE = []
patphoBoost = []
with open(sys.argv[1], 'r') as csvfile:
    test1reader = csv.reader(csvfile, dialect='excel')
    for row in test1reader:
        rows.append(row)
        if row[0] != 'Number of PAT Photons':
            numpatpho.append(int(row[0]))
            patphopt.append(float(row[1]))
            patphoM.append(float(row[2]))
            patphoE.append(float(row[3]))
            patphoBoost.append(float(row[4]))
    numrows = test1reader.line_num

# for i in range(numrows-1):
#     H_clusters.Fill(patphopt[i], numpatpho[i])

# C = TCanvas()
# C.cd()
# H_clusters.Draw("hist")
# C.Print("ClustersvspT1.pdf")

# Emin = float(sys.argv[3]) - 10
# Emax = float(sys.argv[4]) + 10

# # ~~ Dot Plotting ~~

# plt.style.use('_mpl-gallery')
# fig, ax = plt.subplots()
# for i in range(numrows-1):
#     if numpatpho[i] <= 2 and sys.argv[3] != '-' and patphopt[i] != 0:
#         ax.plot(patphopt[i]/float(sys.argv[3]), numpatpho[i], 'o-', linewidth=1)
#     elif numpatpho[i] <= 2 and patphopt[i] != 0:
#         ax.plot(patphopt[i], numpatpho[i], 'o-', linewidth=1)
# ax.set(xlim=(Emin,Emax), xticks=np.arange(Emin,Emax,float(sys.argv[5])/10), ylim=(0, 3), yticks=np.arange(1,3))
# # plt.show()
# if sys.argv[3] != '-':
#     ax.set_title('Number of PAT photons vs. Higgs pT/mass')
#     ax.set_xlabel('Higgs pT/mass')
#     ax.set_ylabel('# Photons')
# else:
#     ax.set_title('Number of PAT photons vs. Higgs pT')
#     ax.set_xlabel('Higgs pT (GeV)')
#     ax.set_ylabel('# Photons')

# plt.tight_layout()
# fig.set_size_inches(9, 5, forward=True)
# # plt.show()
# plt.savefig(sys.argv[2], format='pdf')

# ~~ Histogramming ~~

NBinsX = 0
NBinsY = 0
binpT = [0]
binclust = [0]
counterX = float(sys.argv[3])
counterY = 0
Xbinsizes = [11]
Xbinincrements = [float(sys.argv[4])/10]
# Xbinincrements = [500/10]
# Xbinincrements = [0.1/10]
Ybinsizes = [1, 6]
Ybinincrements = [0.5, 1]
# ~~ Defining X and Y bins ~~
for j in range(len(Xbinsizes)):
    for i in range(Xbinsizes[j]):
        counterX += Xbinincrements[j]
        NBinsX += 1
        binpT.append(counterX)
for j in range(len(Ybinsizes)):
    for i in range(Ybinsizes[j]):
        counterY += Ybinincrements[j]
        NBinsY += 1
        binclust.append(counterY)

H_efficiency = TH2F("Efficiency", "2 Photon Pair Mass: 3 Photon Pair Mass", (NBinsX-1), array('d', binpT), (NBinsY-1), array('d', binclust))
# H_efficiencycount = TH2F("Efficiencycount", "2 Photon Pair Mass: 3 Photon Pair Mass", (NBinsX-1), array('d', binpT), (NBinsY-1), array('d', binclust))
# H_smallEpct = TH2F("Small Mass Percent", "2 Photon Pair Mass: 3 Photon Pair Mass", (NBinsX-1), array('d', binpT), (NBinsY-1), array('d', binclust))

# rows = []
# largemasses = []
# smallmasses = []
# n_events = []
# n_3phoevents = []
# smallEpct = []
# efficiency = []
# with open('efficiencyinfo1.csv', 'r') as csvfile:
#     test1reader = csv.reader(csvfile, dialect='excel')
#     for row in test1reader:
#         rows.append(row)
#         if row[0] != 'Large Mass':
#             largemasses.append(float(row[0]))
#             smallmasses.append(row[1])
#             n_events.append(float(row[2]))
#             n_3phoevents.append(float(row[3]))
#             smallEpct.append(float(row[4]))
#             efficiency.append(float(row[5]))
#     numrows = test1reader.line_num
# for i in range(numrows-1):
#     tmpstr1 = smallmasses[i]
#     tmpstr2 = tmpstr1.replace('p','.')
#     smallmasses[i] = float(tmpstr2)

# Next step: average efficiency values over the ranges of the bin widths, then fill.


# print(NBinsX)
# print(binpT)
summingvar = patphoBoost
bincounter = [0] * (NBinsX)
# print(bincounter)
# print(len(bincounter))
# print(summingvar[3], patphopt[3], patphoM[3])


countertot = 0
for i in range(numrows-1):
    if numpatpho[i] <= 5: # and patphopt[i] != 0:
        countertot += 1
        for j in range(NBinsX):
            if binpT[j] <= summingvar[i] < binpT[j+1]:
                bincounter[j] += 1
for i in range(numrows-1):
    if numpatpho[i] <= 5:
        for j in range(NBinsX):
            if binpT[j] <= summingvar[i] < binpT[j+1]:
                H_efficiency.Fill(summingvar[i], numpatpho[i], 1/bincounter[j])
            # # if doing Mass/Energy, need to remove any zeros: (not working ): )
            # if patphoE[i] == 0:
            #     countertot -= 1
            #     bincounter[j] -= (1/NBinsX)
            #     continue
            # if binpT[j] <= patphoM[i]/patphoE[i] < binpT[j+1] and bincounter[j] != 0:
            #     H_efficiency.Fill(patphoM[i]/patphoE[i], numpatpho[i], 1/bincounter[j])
    #     # H_efficiency.Fill(patphopt[i], numpatpho[i], 1/counter1)
    #     # H_efficiency.Fill(patphoM[i], numpatpho[i], 1/counter1)
    #     if patphoE[i] == 0:
    #         countertot -= 1
    #         continue
    #     else:
    #         H_efficiency.Fill(patphoM[i]/patphoE[i], numpatpho[i], 1/countertot)
    # # H_smallEpct.Fill(smallmasses[i], largemasses[i], smallEpct[i])

# print(bincounter)
# print(countertot)

# for j in range(NBinsY):
#     for i in range(NBinsX):
#         temp1 = H_efficiencycount.GetBinContent(i, j)
#         if temp1 > 1:
#             temp2 = H_efficiency.GetBinContent(i, j)
#             temp3 = H_smallEpct.GetBinContent(i, j)
#             H_efficiency.SetBinContent(i,j, temp2/temp1)
#             # H_smallEpct.SetBinContent(i,j,temp3/temp1)

H_efficiency.SetStats(0)
# H_efficiency.SetTitle("PAT Clusters vs. Higgs pT")
# H_efficiency.SetTitle("PAT Clusters vs. Higgs Mass")
# H_efficiency.SetTitle("PAT Clusters vs. Higgs MoE")
H_efficiency.SetTitle("PAT Clusters vs. Higgs Boost")
# H_efficiency.SetTitle("CLUE Clusters vs. Higgs pT")
# H_efficiency.SetTitle("CLUE Clusters vs. Higgs Mass")
# H_efficiency.SetTitle("CLUE Clusters vs. Higgs MoE")
gStyle.SetPalette(109)
gStyle.SetPaintTextFormat("4.2f")
C1 = TCanvas()
C1.cd()
H_efficiency.Draw("hist")
H_efficiency.Draw("colztext")
# C1.SetLogx()
# C1.SetLogy()
xAxis = H_efficiency.GetXaxis()
# xAxis.SetTitle("Higgs pT (GeV)")
# xAxis.SetTitle("Higgs Mass (GeV)")
# xAxis.SetTitle("Higgs MoE")
xAxis.SetTitle("Higgs Boost")
xAxis.CenterTitle(kTRUE)
# gStyle.SetNdivisions(n=len(smallmasses), axis='x')
# gStyle.SetPadTickX(smallmasses)


yAxis = H_efficiency.GetYaxis()
yAxis.SetTitle("PAT Clusters")
# yAxis.SetTitle("CLUE Clusters")
yAxis.CenterTitle(kTRUE)
C1.Print(sys.argv[2])
# C1.Print("efficiencymap7.root")
