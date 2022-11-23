from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'DYToEE_Pythia6_GEN'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'gen_DYToEE_Pythia6.py'
config.Data.outputPrimaryDataset = 'DYToEE_13TeV_Pythia6'

config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 3000
NJOBS = 1000
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
config.Data.publication = True 
config.Site.storageSite = 'T3_US_FNALLPC'
