import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/QCD_Pt-15to7000_TuneCP5_Flat_13TeV_pythia8/QCD_Pt-15to7000_TuneCP5_13TeV_pythia8_Flat_DIGI-RECO/230210_124030/0000/digiToRecoStep_753.root' #QCD
inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/TTToHadronic_TuneCP5_13TeV_powheg-pythia8/TTToHadronic_TuneCP5_13TeV_powheg-pythia8_DIGI-RECO/230225_224742/0000/digiToReco_GluGluHtoTauTau_126.root' #TTbar hadronic
#inputFiles_='file:/eos/uscms/store/group/lpcml/rchudasa/MCGeneration/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8/WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8_DIGI-RECO/230218_132746/0000/digiToRecoStep_318.root' #WjetstoLNu
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia8/GluGluHToTauTau_Hadronic_M125_13TeV_powheg_pythia_DIGI-RECO/230225_224255/0000/digiToReco_GluGluHtoTauTau_54.root'
#inputFiles_='file:/eos/uscms/store/group/lpcml/bbbam/MCGeneration/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_M-50_13TeV-powheg_pythia8_DIGI-RECO/230228_103050/0000/digiToReco_GluGluHtoTauTau_221.root'

maxEvents_=2
#maxEvents_=-1
skipEvents_=0#
#outputFile_='qcd.root'
outputFile_='ttbar.root'
#outputFile_='WJets.root'

cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print '%s'%cmd
os.system(cmd)
