from CRABClient.UserUtilities import config
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters

CFG = 'HToTauTau_Hadronic_m3p6To15_pT0To200_ctau0To3_eta0To1p4_pythia8_noPU_cfg'
#CFG = 'HToTauTau_gun'


# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName>
# To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>

# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea = 'crab_MC'
config.General.requestName = 'HToTauTau_m3p6To15_pT20To200_ctau0To3_eta0To1p4_pythia8_biased'
config.General.transferOutputs = True
config.General.transferLogs = False

# CMS cfg file goes here:
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = '%s.py'%(CFG) # cms cfg file for generating events
config.JobType.maxMemoryMB = 2800

# Define units per job here:
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 200 # units: as defined by config.Data.splitting
#config.Data.totalUnits = 500000 # total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
config.Data.totalUnits = 3000000 # test production
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Site.storageSite = 'T3_US_FNALLPC'
#config.Data.outLFNDirBase = '/store/group/lpcml/' # add your username as subdirectory
config.Data.outLFNDirBase = '/store/user/ddicroce/'
config.Data.outputPrimaryDataset = 'Run2018_GENtoAODSIM_biased'
config.Data.outputDatasetTag = config.General.requestName
