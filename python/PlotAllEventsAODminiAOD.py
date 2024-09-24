# Standard Includes: ----------------------------------------#
import os
import array
from array import *
import glob
import math
import ROOT
from ROOT import *
import sys
import itertools
from itertools import *
from optparse import OptionParser
# Obviously can only be run in a CMSSW Framework
from DataFormats.FWLite import * 
from HLTrigger import *
# -----------------------------------------------------------#

# Plots all events in miniAOD and AOD

# Command line: python3 PlotAllEventsAODminiAOD.py MiniAOD PAT AOD RECO 100 File.pdf
# arguments: MiniAOD file, PAT, AOD file, RECO, # events, graph name

def HardGet(e, L, H): # shorthand def for getting collection from event
	e.getByLabel(L, H)
	if H.isValid() and len(H.product()) > 0: 
		return H.product()
	return False

# Define Histograms
# H_barrelmini = TH2F("minibarrel", "eta:phi", 170, -86., 86., 360, 0., 360.)
# H_barrelaod = TH2F("aodbarrel", "eta:phi", 170, -86., 86., 360, 0., 360.)
H_barrelmini = TH2F("minibarrel", "eta:phi", 195, -97., 97., 374, -16., 378.)
H_barrelaod = TH2F("aodbarrel", "eta:phi", 195, -97., 97., 374, -16., 378.)
# H_barrel = TH2F("barrel", "eta:phi", 172, -86., 86., 360, 0., 360.)

minifile = glob.glob( sys.argv[1]+"*root" ) #will read in the root files, depending on command line argument 1
print(minifile)
minievents = Events(minifile) #format so it can read later

aodfile = glob.glob( sys.argv[3]+"*root" ) #will read in the root files, depending on command line argument 1
print(aodfile)
aodevents = Events(aodfile) #format so it can read later


#labels and handles (tell CMSSW what the collections are called). You'll need to modify this to match what the EventDump says.
HHb = Handle("edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >")
HLbmini = ("reducedEgamma", "reducedEBRecHits", sys.argv[2])
HLbaod = ("reducedEcalRecHitsEB", "", sys.argv[4]) 

n=0
for event in minievents:
    # print(event)
    n+=1
    if not sys.argv[3] == "-1":
        if n > int(sys.argv[5]): break
    
    Hits = HardGet(event, HLbmini, HHb) # here we use my functiono to load the collection
    if Hits == False: continue # if it fails, skip event
    for h in Hits:
        eeDI = EBDetId(h.detid()) # CAST the iterator as a CMSSW barrel hit
        # print(">" + str(h.energy()) + " ("+str(eeDI.ieta())+", "+str(eeDI.iphi())+")")
        # H_barrel.Fill(eeDI.ieta(), eeDI.iphi(), h.energy())
        H_barrelmini.Fill(eeDI.ieta(), eeDI.iphi(), h.energy())

n=0
for event in aodevents:
    # print(event)
    n+=1
    if not sys.argv[6] == "-1":
        if n > int(sys.argv[5]): break
    
    Hits = HardGet(event, HLbaod, HHb) # here we use my functiono to load the collection
    if Hits == False: continue # if it fails, skip event
    for h in Hits:
        eeDI = EBDetId(h.detid()) # CAST the iterator as a CMSSW barrel hit
        # print(">" + str(h.energy()) + " ("+str(eeDI.ieta())+", "+str(eeDI.iphi())+")")
        # H_barrel.Fill(eeDI.ieta(), eeDI.iphi(), h.energy())
        H_barrelaod.Fill(eeDI.ieta(), eeDI.iphi(), h.energy())

C1 = TCanvas()

# H_barrelaod.SetFillColor(kInvertedDarkBodyRadiator)
C1.cd(1)
H_barrelaod.SetStats(0)
H_barrelaod.SetTitle("MiniAOD and AOD ECal Events")
H_barrelaod.Draw("col")
ex1 = TExec("ex1", "gStyle->SetPalette(106)")
ex1.Draw()
gPad.Update()

# H_barrelmini.SetFillColor(kValentine)
C1.cd(2)
H_barrelmini.SetStats(0)
H_barrelmini.SetTitle("MiniAOD and AOD ECal Events")
H_barrelmini.Draw("SAME col")
ex2 = TExec("ex2", "gStyle->SetPalette(109)")
ex2.Draw()
gPad.Update()

xAxis = H_barrelaod.GetXaxis()
xAxis.SetTitle("Eta")
xAxis.CenterTitle(kTRUE)

yAxis = H_barrelaod.GetYaxis()
yAxis.SetTitle("Phi")
yAxis.CenterTitle(kTRUE)

C1.Print(str(sys.argv[6]))
# C1.Print("AllEventsAODminiAODnewgentest1.root")
