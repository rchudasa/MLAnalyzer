from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'DYToEE_Pythia6_SIM'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'sim_DYToEE_Pythia6.py'
#config.Data.outputPrimaryDataset = 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.inputDataset ='/DYToEE_13TeV_Pythia6/lpcml-crab_DYToEE_Pythia6_GEN-6b3b91cdafc28c0ec152bb6300e15779/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
#config.Data.totalUnits = 1500

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
#config.Data.outLFNDirBase = '/store/user/rchudasa/'
config.Data.publication = True 
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
