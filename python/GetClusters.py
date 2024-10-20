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

# Goal: Import clustering and pt info into a csv file

# Command line: python3 python/GetClusters.py MiniAOD PAT 100 1.0 clusters.csv
# arguments: miniaod file minus .root, PAT, # events, mass, name of csv file

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
# photoken = EDGetToken("std::vector<pat::Photon>")



minievents = GetEvents(sys.argv[1])
mass = float(sys.argv[4])

n=0
for event in minievents:
	n+=1
	numpatpho = 0
	numgenpho = 0
	# patphopt = []
	higgspt = 0
	if not sys.argv[3] == "-1":
		if n > int(sys.argv[3]):
			break
	# # event.getByLabel(HLbpho, HHbpho)
	# # print("PAT photons detected: ", (len(HHbpho.size())))
	# event.getByLabel(HLbgen, HHbgen)
	# # print("Gen photons detected: ", len(HHbgen.product()))
	# print(HHbgen.product())
	
    #event.getByLabel(edm::InputTage, edm::Handle<T>)
	#event.getByToken(edm::EDGetTokenT, edm::Handle<T>)
	# gamma = event.getByToken(photoken, HHbpho)
	# print(gamma.size())
	event.getByLabel(HLbpho, HHbpho)
	totvecpho = []
	vecphotemp = TLorentzVector()
	higgsvec = TLorentzVector()
	for i,pho in enumerate(HHbpho.product()):
		# if pho.pt() > 14 and pho.hadTowOverEm() < 0.15:
		# # if pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3:
		# 	continue
		# print(i, pho.pt(), pho.superCluster().eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta())
		# print(pho.photonID('cutBasedPhotonID-Fall17-94X-V2-tight'))
		# if pho.photonID('cutBasedPhotonID-Fall17-94X-V2-loose') == True:
		numpatpho += 1
		# print(pho.pt(), pho.eta(), pho.phi(), pho.energy())
		vecphotemp.SetPtEtaPhiM(pho.pt(), pho.eta(), pho.phi(), pho.mass())
		higgsvec += vecphotemp
		# totvecpho.append(vecphotemp)
		# higgspt += pho.pt()
	# for i in range(numpatpho):
	# 	higgsvec += totvecpho[i]
	higgspt = higgsvec.Pt()
	higgsE = higgsvec.E()
	higgsboostvec = higgsvec.BoostVector()
	higgsboost = higgsboostvec.Mag()
	# print(numpatpho)
	# if numpatpho >= 2:
	# print(higgspt, higgsvec.Eta(), higgsvec.Phi(), higgsvec.E(), higgsvec.M())

	genvecphotemp = TLorentzVector()
	genhiggsvec = TLorentzVector()
	event.getByLabel(HLbgen, HHbgen)
	for i,gen in enumerate(HHbgen.product()):
		if gen.pdgId() == 22:
			numgenpho += 1
			# genvecphotemp.SetPtEtaPhiE(gen.pt(), gen.eta(), gen.phi(), gen.energy())
			# genhiggsvec += genvecphotemp
		if gen.pdgId() == 25:
			# print(gen.mass())
			higgsmass = gen.mass()
			# higgsmoe = higgsmass/higgsE
		# if gen.pdgId() == 35:
		# 	print(gen.mass())
	# print(numgenpho)
	# print(genhiggsvec.Pt(), genhiggsvec.Eta(), genhiggsvec.Phi(), genhiggsvec.E(), genhiggsvec.M())
	# print("Next event")
	
	if numgenpho == 2: # and (mass-(mass/10) < higgsvec.M() < mass+(mass/10)):
		# print("found")
		with open(sys.argv[5], 'a', newline='') as csvfile:
			test2writer = csv.writer(csvfile, dialect='excel')
			test2writer.writerow([numpatpho] + [higgspt] + [higgsmass] + [higgsE] + [higgsboost])

	



