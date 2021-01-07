import os

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_ ='file:/eos/cms/store/user/mandrews/ML/FEVTDEBUG/h24gamma_1j_10K_100MeV_FEVTDEBUG_2016_25ns_Moriond17MC_PoissonOOTPU/180109_112606/0000/step_full_1.root'
#inputFiles_ ='file:/uscms/home/jda102/nobackup/BTaggingML/CMSSW_8_0_30/src/step_AODSIM.root'
#inputFiles_ = 'file:/eos/uscms/store/user/jda102/AODSIM/QCD_Pt_30_70_13TeV_TuneCUETP8M1_noPU_AODSIM/181019_231226/0000/step_AODSIM_83.root'
#inputFiles_='/store/group/lpcml/mandrews/AODSIM/QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU_AODSIM/180809_215549/0000/step_full_1.root'
#inputFiles_='file:../HToTauTau_m3p6To15_pT20To200_ctau0To3_eta0To1p4_pythia8_biased_v1.root'
#inputFiles_='file:../DoubleTauPt5_50_pythia8_noPU.root'
#inputFiles_='file:../ZTauTau_Tauola_All_hadronic_13TeV_TuneCUETP8M1.root'
#inputFiles_='file:../HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_pythia8.root'
#inputFiles_='file:../DoubleTauPt50_200_pythia8_noPU.root'
#inputFiles_='file:../WJetsToLNu_13TeV_pTHat20To50.root'
#inputFiles_='file:../WJetsToLNu_13TeV_pTHat50To200.root'
#inputFiles_='file:../QCD_Pt-50to200_TuneCUETP8M1_13TeV_pythia8.root'
#inputFiles_='file:../Tau_Run2016B-07Aug17_ver1-v1_MINIAOD_12172182-C499-E711-A87A-0026B9533C1C.root'
inputFiles_='file:../Run2018A_Tau_AOD_12Nov2019_UL2018_1E64EE69-D0C6-3E47-9724-FDF606418876.root'
#inputFiles_='file:../HTauTau_Tune4C_13TeV_LHE_pythia8_Tauola.root'
#inputFiles_='file:../TTbar_TuneCUETP8M1_13TeV_pythia8.root'

maxEvents_=-1
skipEvents_=0#


cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
