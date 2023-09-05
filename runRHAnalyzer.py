import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:../HToEleEle_m100To2000_pT20To120_ctau0To3_eta0To1p4_pythia8_noPU_biased.root'
inputFiles_='file:../HTauTau_Hadronic_Tune4C_13TeV_LHE_pythia8_Tauola_PU2018-RunIISummer20UL18MiniAOD.root'

maxEvents_=-1
skipEvents_=0#
processTask_='tau_classification'
processIsDebug_=True
processIsData_=False
processIsSignal_=True

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d processTask=%s processIsData=%s processIsSignal=%s processIsDebug=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,processTask_,processIsData_,processIsSignal_,processIsDebug_)
print '%s'%cmd
os.system(cmd)
