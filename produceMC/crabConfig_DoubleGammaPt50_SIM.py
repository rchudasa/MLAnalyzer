from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'DoubleGammaPt50_AODSIM'
config.General.workArea = 'crab_inference'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'sim_DoubleGammaPt50_Pythia8.py'
#config.Data.outputPrimaryDataset = 'GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.inputDataset ='/DoubleGammaPt50_Pythia8/phys_heavyions-crab_DoubleGammaPt50_GEN-69a5e0f6d1225b384d191f735ec50b21/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
#config.Data.totalUnits = 1500

#config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/e2e'
config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
