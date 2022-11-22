from CRABClient.UserUtilities import config#, getUsernameFromSiteDB
config = config()
# See parameter defintions here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3ConfigurationFile#CRAB_configuration_parameters

#idx = '00000'
#CFG = 'QCD_Pt_80_170_%s'%idx
#CFG = 'HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_biased_v2'
CFG = 'GJet_Pt-20to40'

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
config.JobType.psetName = 'RecHitAnalyzer/python/ConfFile_cfg.py' # analyzer cfg file
#config.JobType.maxMemoryMB = 2800

# Define input and units per job here:
#config.Data.userInputFiles = open('MLAnalyzer/list_production.txt'%idx).readlines()
config.Data.inputDataset ='/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIISummer19UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 # units: as defined by config.Data.splitting
config.Data.totalUnits = -1 # -1: all inputs. total jobs submitted = totalUnits / unitsPerJob. cap of 10k jobs per submission
#config.Data.totalUnits = 10 # test production
config.Data.publication = False

# Output files will be stored in config.Site.storageSite at directory:
# <config.Data.outLFNDirBase>/<config.Data.outputPrimaryDataset>/<config.Data.outputDatasetTag>/
#config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/e2e' # add your username as subdirectory
#config.Data.outputPrimaryDataset = 'HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_biased'
#config.Data.outputPrimaryDataset = 'DYToEE_ntuples'
config.Data.outputDatasetTag = config.General.requestName
