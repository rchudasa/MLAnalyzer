import os

cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/JetHT-AOD/JetHT/JetHT_Run2018A-v1_RAW-AOD-multithread/240621_165042/0000/RAW2DIGI_L1Reco_RECO_trackRechitsAdded_jetHT_56.root'

#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/Tau/Tau_Run2018D_RAW-AOD-multithread/240208_082349/0000/RAW2DIGI_L1Reco_RECO_trackRechitsAdded_2.root'
#inputFiles_='file:../dataReco/RAW2DIGI_L1Reco_RECO_trackRechitsAdded.root'
#inputFiles_='file:../Run2018A_Tau_AOD_12Nov2019_UL2018_1E64EE69-D0C6-3E47-9724-FDF606418876.root'
#inputFiles_='file:../Run2018A_Tau_AOD_12Nov2019_UL2018_D0938794-D6B0-594C-BB05-5CEE65208347.root'

maxEvents_=10
skipEvents_=0#
outputFile_='data_event-AOD_triggerTree.root'

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print '%s'%cmd
os.system(cmd)
