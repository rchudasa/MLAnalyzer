import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:../MCProduction/E2E-TauClassifier/digiToRecoStep.root'
#inputFiles_='file:../MCProduction/E2E-TauClassifier/digiToRecoStep_multiThread.root'

#maxEvents_=200
maxEvents_=10
skipEvents_=0#

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
