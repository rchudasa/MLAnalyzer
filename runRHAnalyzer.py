import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:../HToEleEle_m100To2000_pT20To120_ctau0To3_eta0To1p4_pythia8_noPU_biased.root'
#inputFiles_='file:../HTauTau_Hadronic_Tune4C_13TeV_LHE_pythia8_Tauola_PU2018-RunIISummer20UL18MiniAOD.root'
#inputFiles_='file:/eos/uscms/store/group/lpcml/ddicroce/Samples/HTauTau_13TeV-RunIISummer20UL18MiniAOD_validation/0000/HTauTau_Hadronic_Tune4C_13TeV_LHE_pythia8_Tauola_PU2018-RunIISummer20UL18MiniAOD_99.root'
#inputFiles_='file:/uscms/home/rchudasa/nobackup/miniAOD_checks/CMSSW_10_6_25/src/miniAOD_production/DYToTauTau_MiniAOD.root'
inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8_MINIAODSIM/230908_060048/0000/DYToTauTau_MiniAOD_889.root'

maxEvents_=100
skipEvents_=0#
processTask_='tau_classification'
processIsDebug_=True
processIsData_=False
processIsSignal_=False
processIsW_=True
outputFile_="WJets_myProduction.root"
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d processTask=%s processIsData=%s processIsSignal=%s processIsDebug=%s processIsW=%s outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,processTask_,processIsData_,processIsSignal_,processIsDebug_,processIsW_,outputFile_)
print '%s'%cmd
os.system(cmd)
