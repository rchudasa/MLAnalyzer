from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'QCDEmEnriched_Pt-30to50_SIM-v2'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'sim_GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.py'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.inputDataset ='/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8/phys_heavyions-crab_QCD_Pt-30to50_EMEnriched_GEN-65d7ea1a490a5182af464ee341e868ac/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
#config.Data.outLFNDirBase = '/store/user/rchudasa/'
config.Data.publication = True 
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
