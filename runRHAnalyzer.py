import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:../HToAAToTauTau_Hadronic_M5_13TeV_2018.root'
inputFiles_='file:../HToTauTau_m3p6To17_pT20To200_ctau0To3_eta0To1p4_pythia8_biased_v2.root'
#inputFiles_='file:../HToTauTau_m3p6To15_pT20To200_ctau0To3_eta0To1p4_pythia8_biased_v2.root'

maxEvents_=-1
skipEvents_=0#


cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
