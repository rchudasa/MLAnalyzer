from CRABClient.UserUtilities import config
config = config()

CFG = 'DYToTauTau_ntuples_miniAOD'

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
#config.Data.userInputFiles = open('DYToTauTau.txt').readlines()
#config.Data.inputDataset ='/DYToTauTau_M-50_13TeV-powheg_pythia8/lpcml-DYToTauTau_M-50_13TeV-powheg_pythia8_DIGI-RECO-2941225365836e99c29b2389a87cd4a2/USER'
config.Data.inputDataset ='/DYToTauTau_M-50_13TeV-powheg_pythia8/lpcml-DYToTauTau_M-50_13TeV-powheg_pythia8_MINIAODSIM-7ee272c68135d2b13931732a52559311/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2 # units: as defined by config.Data.splitting
config.Data.totalUnits = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/NTuples' # add your username as subdirectory
#config.Data.outputPrimaryDataset = 'DYToTauTau_M-50_13TeV-powheg_pythia8'#if providing text file input then primary dataset name is required
config.Data.outputDatasetTag = config.General.requestName
