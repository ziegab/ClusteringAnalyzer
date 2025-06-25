import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018
import sys

process = cms.Process("Clus")

SkipEvent = cms.untracked.vstring('ProductNotFound')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'

nevents = 500
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents) )

process.source = cms.Source("PoolSource",
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events0p1Ma0p01_101p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events0p2Ma0p01_201p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_500events0p4Ma0p01_401p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events0p6Ma0p01_601p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events0p8Ma0p01_801p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events2p0Ma1p0_2001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events3p0Ma1p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events4p0Ma1p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events5p0Ma1p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events6p0Ma1p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events8p0Ma10p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events10p0Ma10p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events15p0Ma10p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events20p0Ma10p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamples/AtoGG_1000events25p0Ma10p0_3001p0EPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/DiphotonSamplesNoResonance/AtoGG_1000events1p0Ma0p01_1001p0EPileupNoResv1_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/eos/user/g/gziemyte/AtoGG_5000events0p01_0p05MoE10p0_200p0EPileupNoResv55_MiniAOD.root')
                                fileNames = cms.untracked.vstring('file:'+sys.argv[2])
                            )
# Input histogram file to be filled maybe
process.TFileService = cms.Service("TFileService", 
                                  #  fileName = cms.string('hist_ClusteringAnalyzer_AtoGG_500events1.0Ma1000.10000pTPU.root')
                                  #  fileName = cms.string('hist_AtoGG_1000events25p0Ma2dc15rhoc3delt.root')
                                  fileName = cms.string('hist_AtoGG_'+str(nevents)+'events_'+sys.argv[4]+'_MoE2dc15rhoc3deltnoresv'+sys.argv[3]+'.root')
                                  #  fileName = cms.string('histplots_AtoGG_100events0p1Ma2dc15rhoc3delt.root')
                                  #  fileName = cms.string('hist_ClusteringAnalyzer_AtoGG500events10.0Ma100.3000pTPU.root')
                                  #  fileName = cms.string('hist_ClusteringAnalyzer_1Dinfo')
                                #    fileName = cms.string('hist_ClusteringAnalyzer_BkkToGRadionToGGG_M1-500_R0-5.root')
)
eventN = list(range(1, (nevents+1)))
process.clus = cms.EDAnalyzer('ClusteringAnalyzer',
   ecalRechitsEB = cms.InputTag('reducedEgamma','reducedEBRecHits'),
   ecalRechitsEE = cms.InputTag('reducedEgamma','reducedEERecHits'),
   ecalRechitsES = cms.InputTag('reducedEgamma','reducedESRecHits'),
   genParticles = cms.InputTag('packedGenParticles'),
   genRecoParticles = cms.InputTag('prunedGenParticles'),
   photons = cms.InputTag('slimmedPhotons'),
  #  electrons = cms.InputTag('slimmedElectrons'),
   EventsToScan = cms.vint32(eventN),
                              )

process.p = cms.Path(process.clus)
