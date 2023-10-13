from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'DYToTauTau_ntuples-AOD'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runRHAnalyzer_DYToTauTau.py'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
config.Data.userInputFiles = open('DYToTauTau_M-50_13TeV-powheg_pythia8_AOD.txt').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2 
config.Data.outputPrimaryDataset = 'DYToTauTau_M-50_13TeV-powheg_pythia8'

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/NTuples'
#config.Data.publication = True 
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
