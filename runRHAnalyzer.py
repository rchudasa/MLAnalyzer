import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:../../../CMSSW_10_2_20_UL/src/HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_pythia8.root'
#inputFiles_='file:../../../CMSSW_10_2_20_UL/src/HTauTau_Tune4C_13TeV_LHE_pythia8_Tauola.root'
#inputFiles_='file:../../../CMSSW_10_2_20_UL/src/ZTauTau_Tauola_All_hadronic_13TeV_TuneCUETP8M1.root'

maxEvents_=-1
skipEvents_=0#


cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
