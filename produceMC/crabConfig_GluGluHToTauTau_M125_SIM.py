from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'GluGluHToTauTau_AODSIM'
config.General.workArea = 'crab_inference'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'sim_GluGluHToTauTau_M125_13TeV_powheg_pythia8.py'
#config.Data.outputPrimaryDataset = 'GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.inputDataset ='/GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8/phys_heavyions-crab_GluGluHToTauTau_M125_GEN-LHEoutput-5cf27f8dc192866c7528a9c0ed7a4cee/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
#config.Data.totalUnits = 1500

#config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/e2e'
config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
