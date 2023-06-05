import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/eventGenerationChecks/digiToRecoStep_753_QCD_Pt-15to7000.root' #QCD
#inputFiles_='file:/eos/cms/store/group/phys_heavyions/rchudasa/e2e/eventGenerationChecks/digiToRecoStep_318_Wjets.root' #QCD
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/TTToHadronic_TuneCP5_13TeV_powheg-pythia8/TTToHadronic_TuneCP5_13TeV_powheg-pythia8_DIGI-RECO/230225_224742/0000/digiToReco_GluGluHtoTauTau_126.root' #TTbar hadronic
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8/GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia_DIGI-RECO/230225_224255/0000/digiToReco_GluGluHtoTauTau_54.root'
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_M-50_13TeV-powheg_pythia8_DIGI-RECO/230228_103050/0000/digiToReco_GluGluHtoTauTau_221.root'#DY To TauTau
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8_DIGI-RECO/230218_132746/0000/digiToRecoStep_426.root'#WJets
inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/DYToEE_M-50_13TeV-powheg_pythia8/DYToEE_M-50_13TeV-powheg_pythia8_DIGI-RECO/230211_163838/0000/digiToReco_withPileup_621.root'#DYToee
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/QCD_Pt-30to50_EMEnriched_TuneCP5_13TeV_pythia8/QCD_Pt-30to50_EMEnriched_DIGI-RECO/230210_122831/0000/digiToReco_withPileup_908.root'#QCDEmEnriched

maxEvents_=20
#maxEvents_=-1
skipEvents_=0#
#outputFile_='qcd.root'
#outputFile_='ttbar_tauCode.root'
#outputFile_='WJets.root'
outputFile_='dyToEE.root'
#outputFile_='acd_EmEnriched.root'

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print '%s'%cmd
os.system(cmd)
