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
import uproot
# -----------------------------------------------------------#

file_dirdc2rhoc15delt3 = '/afs/cern.ch/user/g/gziemyte/private/CMSSW_13_2_4/src/test/clusteringanalyzer/dc2rhoc15delt20'
print(file_dirdc2rhoc15delt3)
root_filesdc2rhoc15delt3 = glob.glob(f"{file_dirdc2rhoc15delt3}/*.root")

# ~~ Histogramming ~~

NBinsX = 0
NBinsY = 0
binpT = [0]
binclust = [0]
counterX = 0
counterY = 0
Xbinsizes = [12]
Xbinincrements = [1100/11]
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

H_efficiency = TH2F("PAT Photon Boost Efficiency", "Boost: Fraction of Events", (NBinsX-1), array('d', binpT), (NBinsY-1), array('d', binclust))

totalEBClustsdc2rhoc15delt3 = 0
totalnevents = 0
clustpereventdc2rhoc15delt3counter = 0
fileswithvalidclust = 0
totalenergyratio = 0
for root_file in root_filesdc2rhoc15delt3:
    counterEBClusts = 0
    neventsinboostrange = 0
    countereventswithEBclust = 0
    with uproot.open(root_file) as file:
        tree = file["clus"]
        Events = tree["Events;1"]
        gammavalues = Events["gammaval"].array(library='np')
        genenergy = Events["EH"].array(library='np')
        nPATpho = Events["nphoton"].array(library='np')
        nclust = Events["nClusters"].array(library='np')
        EBClustE = Events["EB_cluster_E"].array(library='np')
        nevents = len(gammavalues)
        totalenergyratiopermasspoint = 0
        for i in range(nevents):
            tempEBclustE = 0
            EBClustEcell = EBClustE[i]
            if 9 < gammavalues[i] < 11:
                neventsinboostrange += 1
                if len(EBClustEcell) > 0:
                    countereventswithEBclust += 1
                    for j in range(len(EBClustEcell)):
                        counterEBClusts += (1)
                        tempEBclustE += EBClustEcell[j]
            EBclustEperevent = tempEBclustE
            EBclustEbyEHiggsperevent = EBclustEperevent/genenergy[i] # energy ratio per event
            totalenergyratiopermasspoint += EBclustEbyEHiggsperevent
    totalnevents += nevents
    # print(counterEBClusts)
    if countereventswithEBclust > 0:
        fileswithvalidclust += 1
        totalenergyratio += (totalenergyratiopermasspoint/countereventswithEBclust)#/len(root_filesdc2rhoc15delt3)
    if neventsinboostrange > 0:
        totalEBClustsdc2rhoc15delt3 += counterEBClusts#/neventsinboostrange
    # print(totalEBClustsdc2rhoc15delt3)
    # clustpereventdc2rhoc15delt3counter += (totalEBClustsdc2rhoc15delt3/nevents)/len(root_filesdc2rhoc15delt3)
totclustdc2rhoc15delt3 = totalEBClustsdc2rhoc15delt3#/fileswithvalidclust#len(root_filesdc2rhoc15delt3)
# print(len(root_filesdc2rhoc15delt3))
# print(totalEBClustsdc2rhoc15delt3)
print(totclustdc2rhoc15delt3)
print(totalenergyratio/fileswithvalidclust)
print(totalnevents)
        # print(EBClustE)
    # print(counter)
    # print(nevents)
    # F = TFile(root_file)
    # T = F.Get("clus/Events;1")
    # for e in T:
    #     gammavals = e.gammaval
    #     nPATpho = e.nphoton
    #     nClust = e.nClusters
    #     EHiggs = e.EH
    #     EClust = e.cluster_E

    #     if EBClustE:
    #         counter += 1
    # print(counter)

