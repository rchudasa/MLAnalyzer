from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters

#idx = '00000'
CFG = 'QCD_Pt-15to7000_ntuples-AOD'

# To submit to crab:
# crab submit -c crabConfig_data.py
# To check job status:
# crab status -d <config.General.workArea>/<config.General.requestName># To resubmit jobs:
# crab resubmit -d <config.General.workArea>/<config.General.requestName>

# Local job directory will be created in:
# <config.General.workArea>/<config.General.requestName>
config.General.workArea = 'crab_projects'
config.General.requestName = CFG
config.General.transferOutputs = True
config.General.transferLogs = True

# CMS cfg file goes here:
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'runRHAnalyzer_QCD_Pt-15to7000.py' # analyzer cfg file
config.JobType.maxMemoryMB = 4000 

# Define input and units per job here:
config.Data.inputDBS = 'phys03'
config.Data.inputDataset ='/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/lpcml-QCD_Pt-15to7000_TuneCP5_13TeV_pythia8_Flat_DIGI-RECO-fc26027cf5f8fe837f33c898e4f0f7f7/USER'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 # units: as defined by config.Data.splitting
config.Data.totalUnits = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
#config.Data.totalUnits = 10 # test production
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/NTuples' # add your username as subdirectory
#config.Data.outputPrimaryDataset = 'QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8'
config.Data.outputDatasetTag = config.General.requestName
