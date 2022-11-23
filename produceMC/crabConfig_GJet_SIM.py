from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'GJet_Pt-20to40_SIM-v3'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'sim_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.py'
#config.Data.outputPrimaryDataset = 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8'

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.inputDataset ='/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/phys_heavyions-crab_GJet_Pt-20to40_GENv3-99d9aa36932aa71f9045d950d2e50dc5/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
config.Data.totalUnits = 1000

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
#config.Data.outLFNDirBase = '/store/user/rchudasa/'
config.Data.publication = True 
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
