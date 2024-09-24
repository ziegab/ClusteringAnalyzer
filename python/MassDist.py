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
# -----------------------------------------------------------#

# Goal: Get mass distributions from gen and pat info

# Command line: python3 python/GetClusters.py MiniAOD PAT 100 hist.pdf
# arguments: miniaod file minus .root, PAT, # events, name of pdf file

def GetEvents(inputfile): # inputfile cannot have .root ending to it (just )
	file = glob.glob( inputfile+"*root" ) #will read in the root files, depending on command line argument 1
	events = Events(file) #format so it can read later
	return events

#labels and handles (tell CMSSW what the collections are called). You'll need to modify this to match what the EventDump says.
HHb = Handle("edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> >")
HLbmini = ("reducedEgamma", "reducedEBRecHits", sys.argv[2])
# HLbaod = ("reducedEcalRecHitsEB", "", sys.argv[5])
HHbgen = Handle("vector<reco::GenParticle>")
HLbgen = ("prunedGenParticles", "", "PAT")
HHbpho = Handle("std::vector<pat::Photon>")
HLbpho = ("slimmedPhotons", "", "PAT")

H_genmass = TH1F("genmass", "mass:events", 100, -20., 500.)
H_patmass = TH1F("patmass", "mass:events", 100, -20., 500.)

minievents = GetEvents(sys.argv[1])

n=0
for event in minievents:
	n+=1
	numpatpho = 0
	numgenpho = 0
	higgspt = 0
	if not sys.argv[3] == "-1":
		if n > int(sys.argv[3]):
			break
	event.getByLabel(HLbpho, HHbpho)
	totvecpho = []
	vecphotemp = TLorentzVector()
	higgsvec = TLorentzVector()
	for i,pho in enumerate(HHbpho.product()):
		numpatpho += 1
		vecphotemp.SetPtEtaPhiM(pho.pt(), pho.eta(), pho.phi(), pho.mass())
		higgsvec += vecphotemp
	higgsmass = higgsvec.M()
	H_patmass.Fill(higgsmass)

	genvecphotemp = TLorentzVector()
	genhiggsvec = TLorentzVector()
	event.getByLabel(HLbgen, HHbgen)
	for i,gen in enumerate(HHbgen.product()):
		if gen.pdgId() == 22:
			numgenpho += 1
		if gen.pdgId() == 25:
			# print(gen.mass())
			H_genmass.Fill(gen.mass())
		# if gen.pdgId() == 35:
		# 	print(gen.mass())

C1 = TCanvas()
C1.cd()
# gStyle.SetOptTitle(kFALSE)
gStyle.SetOptStat(0)
# hs1.Draw("hist")
H_genmass.Draw("PLC PMC hist")
H_patmass.Draw("SAME PLC PMC hist")

legend1 = TLegend(0.9,0.8,0.7,0.9)
legend1.SetHeader("Legend","C")
legend1.AddEntry(H_genmass, "GenPho Higgs Mass")
legend1.AddEntry(H_patmass, "PatPho Higgs Mass")
legend1.Draw()

C1.Print(str(sys.argv[4]))
# C1.Print("AllEventsAODminiAODnewgentest1.root")
