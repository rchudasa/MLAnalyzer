import os

cfg='RecHitAnalyzer/python/SCRegressor_cfg.py'
#cfg='RecHitAnalyzer/python/SCRegressor_cfg_data.py'
#inputFiles_='file:../step_full.root'
inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/RunIISummer20UL18MiniAODv2-106X_upgrade2018_realistic_v16_L1v1-v4/MINIAODSIM/4ABBA1E5-92FE-0045-ADC3-CF7FCC59C969.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/FEVTDEBUG/h24gamma_1j_1M_100MeV_noPU_FEVTDEBUG/180109_233954/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m000/180413_215734/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m100/180413_215901/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m000/180416_105657/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/DoublePi0Pt30To90_pythia8_m000_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180604_160025/0004/step_full_4917.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/DoublePi0Pt50To60_m0To1600_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180620_123051/0000/step_full_filtered_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/DoublePi0Pt50To60_m000_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180622_100748/0000/step_full_1.root'
#inputFiles_='file:/eos/uscms/store/user/mba2012/AODSIM/SinglePi0Pt60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM_m600/180611_092217/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePi0Pt30To50_m0To1600_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180917_151524/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePi0Pt15To100_m0To1600_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/181212_004828/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePi0Pt15To100_m0To1600_pythia8_noPU_AODSIM/190116_160039/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePhotonPt15To100_pythia8_noPU_AODSIM/190219_231402/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePi0Pt15To100_m000_pythia8_noPU_AODSIM/190220_040632/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePi0Pt15To100_m0To1600_pythia8_noPU_AODSIM_mlog_ptexp/190217_185619/0000/step_full_1.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/100000/72CD8425-3B87-E811-AEA5-24BE05CEADA1.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/BC2864DE-5B42-E811-9C51-0025905A6138.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/DYToEE_M-50_NNPDF31_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/544FA228-5845-E811-A88A-782BCB539B14.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/288ED03F-6C42-E811-96EC-0CC47A4C8E98.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/QCD_Pt-30toInf_DoubleEMEnriched_MGG-40to80_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FEE7EA67-5942-E811-9C25-0CC47A4D75EE.root'
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/04B84384-FF41-E811-ACA3-7845C4FC3B18.root'
#inputFiles_='/store/data/Run2017B/DoubleEG/MINIAOD/31Mar2018-v1/80000/FAC7DC8A-3737-E811-8BA7-6CC2173DC380.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/h24gamma_1j_1M_100MeV_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180903_152402/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/h24gamma_1j_1M_400MeV_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180908_013040/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/h24gamma_1j_1M_1GeV_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180902_213131/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/h24gamma_1j_1M_100MeV_PU2017_MINIAODSIM/190409_144816/0000/step_miniaodsim_10.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/h24gamma_1j_1M_200MeV_PU2017_MINIAODSIM/190410_034928/0000/step_miniaodsim_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/h24gamma_1j_1M_400MeV_PU2017_MINIAODSIM/190410_035004/0000/step_miniaodsim_1.root'
#inputFiles_='/store/mc/RunIIFall17DRPremix/SUSYGluGluToHToAA_AToGG_M-1_TuneCP5_13TeV_pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/7E53D72A-AD71-E811-B026-0025902EECDC.root'
#inputFiles_='/store/mc/RunIISummer16MiniAODv2/SUSYGluGluToHToAA_AToGG_M-1_TuneCUETP8M1_13TeV_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/D6C6297F-ECCB-E611-AB1F-D067E5F91DA6.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/h24gamma_1j_1M_400MeV_PU2017_AODSIM/190404_040951/0000/step_aodsim_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/SM2gamma_1j_1M_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180904_071354/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/h22gammaSM_1j_1M_2016_25ns_Moriond17MC_PoissonOOTPU_AODSIM/180831_225854/0000/step_full_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/h24gamma_1j_1M_1GeV_PU2017_MINIAODSIM_ext3/190525_045306/0000/step_miniaodsim_ext3_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/h24gamma_1j_1M_400MeV_PU2017_MINIAODSIM_ext3/190525_045433/0000/step_miniaodsim_ext3_166.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePi0Pt10To100_m0To1600_pythia8_ReAOD_PU2017_MINIAODSIM/190610_033632/0000/step_MINIAODSIMext2_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePi0Pt15To100_m0To1600_pythia8_PU2017_MINIAODSIM/190510_220104/0000/step_MINIAODSIM_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePhotonPt15To100_pythia8_PU2017_MINIAODSIM/190430_014319/0000/step_MINIAODSIM_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoubleElectronPt15To100_pythia8_PU2017_MINIAODSIM/190506_230453/0000/step_MINIAODSIM_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePi0Pt15To100_m0To1600_pythia8_PU2017_MINIAODSIM/190510_220104/0000/step_MINIAODSIM_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePi0Pt15To100_m0To1600_pythia8_PU2017_AODSIM/190426_023521/0000/step_AODSIM_1.root'
#inputFiles_='/store/user/lpcml/mandrews/AODSIM/DoublePi0Pt10To100_m0To1600_pythia8_PU2017_AODSIM/190516_143234/0000/step_AODSIM_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePi0Pt10To100_m0To1600_pythia8_PU2017_MINIAODSIM_ext/190517_125544/0000/step_MINIAODSIMext_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePi0Pt10To100_m0To1600_pythia8_PU2017_MINIAODSIM_ext2/190521_221811/0000/step_MINIAODSIMext2_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePhotonPt15To100_pythia8_PU2017_MINIAODSIM_ext2/190522_041329/0000/step_MINIAODSIMext2_1.root'
#inputFiles_='/store/user/lpcml/mandrews/MINIAODSIM/DoublePhotonPt10To100_pythia8_ReAOD_PU2017_MINIAODSIM/190611_010300/0000/step_MINIAODSIMext2_1.root'
#inputFiles_='/store/mc/RunIISummer16MiniAODv3/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/40000/7629BBE2-6B25-E911-918C-5065F37DD4B2.root'
#inputFiles_='/store//mc/RunIISummer16MiniAODv3/GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2/270000/0A87CC73-BC63-E911-B47C-AC1F6B8AC09E.root'
#inputFiles_='root://cmseos.fnal.gov//store/user/lpcml/mandrews/MINIAODSIM/h24gamma_2016/h24gamma_1j_1M_200MeV_pythia8_PU2016/190919_123113/0000/step_MINIAOD_1.root'
#inputFiles_='/store/mc/RunIISummer16MiniAODv3/GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCUETP8M1_13TeV_Pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v2/90000/82035724-E425-E911-9E1D-1866DA879ED8.root'
#inputFiles_='file:/uscms/home/mba2012/nobackup/GUN/CMSSW_9_4_13/src/step_MINIAODSIMext2.root'
#inputFiles_='file:step_full_filtered.root'
#inputFiles_='file:QCD_Pt-40toInf.root'
#inputFiles_='file:diphoton.root'
#inputFiles_='file:doubleEG.root'
#inputFiles_='file:GJet_Pt20-40.root'
#inputFiles_='file:doubleMuon.root'
#inputFiles_='file:jetHT.root'
#inputFiles_='root://cmseos.fnal.gov/%s'%inputFiles_
#inputFiles_='file:/eos/uscms%s'%inputFiles_
#inputFiles_='/store/mc/RunIIFall17MiniAODv2/QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/10000/0C57DB95-FA41-E811-B58D-008CFAFBEBF2.root'
#root://cmsxrootd-site.fnal.gov

#maxEvents_=-1
maxEvents_=1000
#maxEvents_=10
skipEvents_=0
#outputFile_='output.root'
#inputTag=inputFiles_.strip('file:').strip('_FEVTDEBUG.root')
#inputTag='TEST'

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
#for ievt in range(1):
#if not os.path.isdir(inputTag):
#    os.system('mkdir %s'%(inputTag))
#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,ievt)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
print '%s'%cmd
os.system(cmd)
#os.system('mv c*.eps %s/'%(inputTag))
#    os.system('mv cEB*_%d.eps %s/'%(ievt+1,inputTag))

#os.system('scram b -j8')
