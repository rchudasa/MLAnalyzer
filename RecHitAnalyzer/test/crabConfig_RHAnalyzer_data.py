from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'Tau_Run2018D_ntuples_fromAOD'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runRHAnalyzer_data.py'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.Data.inputDataset ='/Tau/phys_heavyions-Tau_Run2018D_RAW-AOD-multithread-6d72431937b9b31ce737dc74b7b8c511/USER'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
config.Data.totalUnits  = -1 
#config.Data.outputPrimaryDataset = 'DYToTauTau_M-50_13TeV-powheg_pythia8'

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/NTuples'
config.Data.publication = True 
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
