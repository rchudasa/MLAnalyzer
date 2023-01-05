from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8_AODSIM-v2'
config.General.workArea = 'crab_inference'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'sim_TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8.py'
#config.Data.outputPrimaryDataset = 'GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
#config.Data.inputDataset ='/TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8/phys_heavyions-crab_TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8_GEN-LHEoutput-65ea0015a84c80fc9611393edb6611b0/USER'
#config.Data.inputDataset ='/TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8/phys_heavyions-crab_TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8_GEN-LHEoutput-65ea0015a84c80fc9611393edb6611b0/USER'
config.Data.userInputFiles = open('input_TTHadronic_GEN.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
#config.Data.totalUnits = 1500
config.Data.outputPrimaryDataset = 'TTToHadronic_TuneCP5CR1_13TeV-powheg-pythia8_AODSIM-v2'

#config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/e2e'
config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
