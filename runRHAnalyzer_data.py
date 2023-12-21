import os

cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
inputFiles_='file:/uscms/home/rchudasa/nobackup/CMSSW_10_6_4_patch1/src/RAW2DIGI_L1Reco_RECO_trackRechitsAdded.root'
#inputFiles_='file:../dataReco/RAW2DIGI_L1Reco_RECO_trackRechitsAdded.root'
#inputFiles_='file:../Run2018A_Tau_AOD_12Nov2019_UL2018_1E64EE69-D0C6-3E47-9724-FDF606418876.root'
#inputFiles_='file:../Run2018A_Tau_AOD_12Nov2019_UL2018_D0938794-D6B0-594C-BB05-5CEE65208347.root'

maxEvents_=-1
skipEvents_=0#


cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
