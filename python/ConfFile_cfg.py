import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process("Clus")

SkipEvent = cms.untracked.vstring('ProductNotFound')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/g/gagarwal/work/diphoton/clustering/CMSSW_13_2_4/src/test/ClusteringAnalyzer/BkkToGRadionToGGG_M1-500_R0-5_2018_MiniAOD_0.root')
                                # fileNames = cms.untracked.vstring(
                                #     'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_0E988FEC-F3B8-E811-9E5C-0242AC1C0501_55k.root'#,
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_60EAD587-54B9-E811-8206-0242AC1C0505_61k.root',
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_6E101B89-BC68-E811-B965-782BCB38D552_55k.root',
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_A0D1E8AF-4471-E811-B5ED-141877410ACD_300.root',
                                #     # 'root://cmseos.fnal.gov//store/user/cmsdas/2023/short_exercises/jets/ttbar2017/TTJets_2017_D2965D3E-2168-E811-9250-782BCB558E3C_59k.root'
                                # )
                            )
process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string('hist_ClusteringAnalyzer_BkkToGRadionToGGG_M1-500_R0-5.root')
)
eventN = list(range(1, 11))
process.clus = cms.EDAnalyzer('ClusteringAnalyzer',
   ecalRechitsEB = cms.InputTag('reducedEgamma','reducedEBRecHits'),
   ecalRechitsEE = cms.InputTag('reducedEgamma','reducedEERecHits'),
   ecalRechitsES = cms.InputTag('reducedEgamma','reducedESRecHits'),
   genParticles = cms.InputTag('packedGenParticles'),
   EventsToScan = cms.vint32(eventN),
                              )

process.p = cms.Path(process.clus)
