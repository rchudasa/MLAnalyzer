from CRABClient.UserUtilities import config
config = config()
#Process = '3p7'
#Process = '8'
Process = '12'
#Process = '10'

inputProcess_ ={
        #'3p7': "HToAATo4Tau_M_3p7_AOD_withTrigger.txt"
        #'8': "HToAATo4Tau_M_8_AOD_withTrigger.txt",
        '12': "HToAATo4Tau_M_12_AOD_withTrigger.txt",
        #'WJets': "WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8_GEN-SIM.txt",
        #'TTbar': "TTToHadronic_TuneCP5_13TeV_powheg-pythia8_GEN-SIM.txt"
        }.get(Process, None)

outputDataset_ = {
        #'3p7':'HToAATo4Tau_M_3p7_pythia8_2018UL_AOD'
        #'8':'HToAATo4Tau_M_8_pythia8_2018UL_AOD',
        '12':'HToAATo4Tau_M_12_pythia8_2018UL_AOD',
        #'WJets':'WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8',
        #'TTbar':'TTToHadronic_TuneCP5_13TeV_powheg-pythia8'
        }.get(Process, None)

#config.section_('General')
config.General.requestName = '%s_RHEvent-Trigger_ntuples'%Process
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RecHitAnalyzer/python/ConfFile_cfg.py'
config.JobType.maxMemoryMB = 4000

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
config.Data.userInputFiles = open('%s'%inputProcess_).readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 
config.Data.outputPrimaryDataset = outputDataset_ 

config.Data.outLFNDirBase = '/store/group/phys_diffraction/rchudasa/MCGeneration/withTrigger'
config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
config.Data.outputDatasetTag = config.General.requestName
