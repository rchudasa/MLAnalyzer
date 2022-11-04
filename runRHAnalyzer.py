import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/Run2018_GENtoAODSIM/DYToEE/221028_061311/0000/step_GEN2AODSIM_356.root'
#inputFiles_='file:../HToEleEle_m100To2000_pT20To120_ctau0To3_eta0To1p4_pythia8_noPU_biased.root'
#inputFiles_='file:../HToEleEle_m100To8000_pT20To120_ctau0To3_eta0To1p4_pythia8_noPU_biased.root'
#inputFiles_='file:prepareDataset/DoubleTauPt5_50_pythia8_noPU.root'
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/Run2018_GENtoAODSIM/HToAAToTauTau_13TeV_Hadronic_M10/220928_123044/0000/HToAAToTauTau_Hadronic_M10_13TeV_2018_102.root'
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/Run2018_GENtoAODSIM/aToTauTau_m3p6To15_pT30To200_ctau0To3_eta0To1p4_pythia8_noPU/221011_062903/0000/aToTauTau_m3p6To15_pT30To200_ctau0To3_eta0To1p4_pythia8_noPU_56.root'
#inputFiles_='file:../step_GEN2AODSIM.root'
#inputFiles_='root://eospublic.cern.ch//eos/opendata/cms/MonteCarlo2012/Summer12_DR53X/DYToEE_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6/AODSIM/NewG4Phys_PU_RD1_START53_V7N-v1/00000/007BAF65-77C7-E311-AF8F-001E673984FD.root'

#maxEvents_=10
maxEvents_=-1
skipEvents_=0#

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
