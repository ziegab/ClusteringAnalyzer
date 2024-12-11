import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process("Clus")

SkipEvent = cms.untracked.vstring('ProductNotFound')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'

nevents = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nevents) )

process.source = cms.Source("PoolSource",
                                # fileNames = cms.untracked.vstring('file:/afs/crc.nd.edu/user/g/gziemyt2/data_analysis/CMSSW_13_2_4/src/test/clusteringanalyzer/MiniAOD20161.root')
                                # fileNames = cms.untracked.vstring('file:/afs/crc.nd.edu/user/g/gziemyt2/data_analysis/CMSSW_13_2_4/src/test/clusteringanalyzer/tempstore/AtoGG_250events1.0Ma10.02000.0pTPileup_MiniAODv2.root')
                                # fileNames = cms.untracked.vstring('file:/afs/crc.nd.edu/user/g/gziemyt2/data_analysis/CMSSW_13_2_4/src/test/clusteringanalyzer/tempstore/AtoGG_500events1p5Ma1p0_1501p0EPileup_MiniAODvalidationSample.root')
                                fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/g/gziemyte/diphogen/tempstore/AtoGG_1000events0p1Ma0p01_101p0EPileup_MiniAOD.root')
                                  # 'file:/afs/crc.nd.edu/user/g/gziemyt2/data_analysis/CMSSW_13_2_4/src/test/clusteringanalyzer/tempstore/AtoGG_500events1.0Ma10.03000.0pTPileup_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov///store/mc/BuToKJpsi_Toee_Mufilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/RunIIAutumn18MiniAOD-PUPoissonAve20_BParking_102X_upgrade2018_realistic_v15-v2/MINIAODSIM')
                              #   fileNames = cms.untracked.vstring('file:/afs/crc.nd.edu/user/g/gziemyt2/data_analysis/CMSSW_13_2_4/src/test/clusteringanalyzer/rootfiles/AtoGammaGammaFlatMoE_10events_MiniAOD.root')
                                # fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/g/gagarwal/work/diphoton/clustering/CMSSW_13_2_4/src/test/ClusteringAnalyzer/BkkToGRadionToGGG_M1-500_R0-5_2018_MiniAOD_0.root')
                                # fileNames = cms.untracked.vstring(
                                #     'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_0E988FEC-F3B8-E811-9E5C-0242AC1C0501_55k.root'#,
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_60EAD587-54B9-E811-8206-0242AC1C0505_61k.root',
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_6E101B89-BC68-E811-B965-782BCB38D552_55k.root',
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_A0D1E8AF-4471-E811-B5ED-141877410ACD_300.root',
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_D2965D3E-2168-E811-9250-782BCB558E3C_59k.root'
                                # )
                            )
# Input histogram file to be filled maybe
process.TFileService = cms.Service("TFileService", 
                                  #  fileName = cms.string('hist_ClusteringAnalyzer_AtoGG_500events1.0Ma1000.10000pTPU.root')
                                   fileName = cms.string('hist_ClusteringAnalyzer_AtoGG_1000eventslpc0p1Ma.root')
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
