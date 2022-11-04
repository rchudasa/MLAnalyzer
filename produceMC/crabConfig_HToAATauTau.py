from CRABClient.UserUtilities import config
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters

CFG = 'HToAAToTauTau_13TeV_cfg'


# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName>
# To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>

# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea = 'crab_MC'
config.General.requestName = 'HToAAToTauTau_13TeV'
config.General.transferOutputs = True
config.General.transferLogs = False

# CMS cfg file goes here:
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = '%s.py'%(CFG) # cms cfg file for generating events
config.JobType.maxMemoryMB = 2800

# Define units per job here:
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 500 # units: as defined by config.Data.splitting
#config.Data.totalUnits = 500000 # total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
NJOBS = 250
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS
config.Data.publication = True

# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Site.storageSite = 'T2_CH_CERN'
#config.Data.outLFNDirBase = '/store/group/lpcml/' # add your username as subdirectory
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/e2e'
config.Data.outputPrimaryDataset = 'Run2018_GENtoAODSIM'
config.Data.outputDatasetTag = config.General.requestName
