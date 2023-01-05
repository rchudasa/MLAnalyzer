from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'GluGluHToTauTau_M125_GEN'
config.General.workArea = 'crab_inference'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'gen_GluGluHToTauTau_M125_13TeV_powheg_pythia8.py'
config.Data.outputPrimaryDataset = 'GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8'

config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
#config.JobType.eventsPerLumi=100
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob =5000
NJOBS = 3
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/e2e'
#config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/'
#config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
#config.Site.storageSite = 'T3_US_FNALLPC'
