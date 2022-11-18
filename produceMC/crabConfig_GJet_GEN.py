from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'GEN_Run2018'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8.py'
config.Data.outputPrimaryDataset = 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8'

config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 5000
NJOBS = 500
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/e2e'
config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
