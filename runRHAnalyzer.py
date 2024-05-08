import os

cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
#cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_M-50_13TeV-powheg_pythia8_DIGI-RECO-v2/230718_043012/0000/digiToReco_GluGluHtoTauTau_100.root' #QCD
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/eventGenerationChecks/digiToRecoStep_753_QCD_Pt-15to7000.root' #QCD
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/eventGenerationChecks/digiToRecoStep_318_Wjets.root' #QCD
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/gen_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/sim_HToAATo4Tau_Hadronic_tauDR0p4_M10_ctau0To3_eta0To2p4_pythia8_2018UL/230512_123635/0000/SIM_HToAAToTauTau_M10_2018UL_withPU_66.root' #H->AA->4Tau
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_pythia8/GJet_Pt-20to40_DoubleEMEnriched_DIGI-RECO-v2/230226_161718/0000/digiToReco_withPileup_605.root' #GJet
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/TTToHadronic_TuneCP5_13TeV_powheg-pythia8/TTToHadronic_TuneCP5_13TeV_powheg-pythia8_DIGI-RECO/230225_224742/0000/digiToReco_GluGluHtoTauTau_126.root' #TTbar hadronic
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8/GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia_DIGI-RECO/230225_224255/0000/digiToReco_GluGluHtoTauTau_54.root'
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/eventGenerationChecks/digiToReco_DYToTauTau_385.root'#DY To TauTau
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8_DIGI-RECO/230218_132746/0000/digiToRecoStep_1.root'#WJets
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/DYToEE_M-50_13TeV-powheg_pythia8/DYToEE_M-50_13TeV-powheg_pythia8_DIGI-RECO/230211_163838/0000/digiToReco_withPileup_621.root'#DYToee
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8/QCD_Pt-30to50_EMEnriched_DIGI-RECO/230210_122831/0000/digiToReco_withPileup_908.root'#QCDEmEnriched
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/eventGenerationChecks/QCD_Pt-30to50_EMEnriched_digiToReco_withPileup_908.rooot'#QCDEmEnriched
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/eventGenerationChecks/GJet_Pt-20to40_DoubleEMEnriched_digiToReco_withPU_617.root'#GJet
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/Tau_Run2018A_UL/C69BBD7C-DC9B-3849-89CB-444AF555A078.root'#Data
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/officialMC_DYJetsToLL_M-50_TuneCP5_13Tev_RECO/61F13245-CF73-9946-8321-B18051BB8659.root'#official MC
#inputFiles_='file:../PhaseI_TTbar_13TeV_NoPu_RECO_newGT.root'#pixel checks
inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/Tau/Tau_Run2018D_RAW-AOD-multithread/240208_082349/0000/RAW2DIGI_L1Reco_RECO_trackRechitsAdded_2.root'#pixel checks

maxEvents_=10
#maxEvents_=20
#maxEvents_=-1
skipEvents_=0#
#outputFile_='MLAnal_PhaseI_TTbar_13TeVu_trackRefitter.root'
#outputFile_='GJet.root'
#outputFile_='ttbar_secVertex.root'
#outputFile_='DYToTauTau_subJet.root'
#outputFile_='WJets_secVertex.root'
#outputFile_='dyToEE.root'
#outputFile_='acd_EmEnriched.root'
outputFile_='data_Tau-AOD_v2.root'

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print '%s'%cmd
os.system(cmd)
