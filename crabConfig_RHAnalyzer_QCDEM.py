from CRABClient.UserUtilities import config
config = config()

CFG = 'QCD_EMEnriched-ntuples'

# Local job directory will be created in:
config.General.workArea = 'crab_projects'
config.General.requestName = CFG
config.General.transferOutputs = True
config.General.transferLogs = True

# CMS cfg file goes here:
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RecHitAnalyzer/python/ConfFile_cfg.py' # analyzer cfg file
config.JobType.maxMemoryMB = 5000

# Define input and units per job here:
config.Data.inputDBS = 'phys03'
#config.Data.userInputFiles = open('DYToEE.txt').readlines()
config.Data.inputDataset ='/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8/lpcml-QCD_Pt-30to50_EMEnriched_DIGI-RECO-fc26027cf5f8fe837f33c898e4f0f7f7/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 # units: as defined by config.Data.splitting
config.Data.totalUnits = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/NTuples' # add your username as subdirectory
#config.Data.outputPrimaryDataset = 'DYToEE_M-50_13TeV-powheg_pythia8'#if providing text file input then primary dataset name is required
config.Data.outputDatasetTag = config.General.requestName
