import os
from glob import glob
import re
import argparse

parser = argparse.ArgumentParser(description='Run RHAnalyzer')
args = parser.parse_args()

eosDir='/eos/uscms/store/user/ddicroce/Run2018_GENtoAODSIM_v1/DiPi0ToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_pythia8_noPU_v2_cfg/201020_150516'
#eosDir='/eos/cms/store/user/mandrews'
xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
output='HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4'
outdir='root://eoscms.cern.ch/eos/user/d/ddicroce/ML/TauClassifier'

cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_ = ['file:%s'%path for path in glob('%s/FEVTDEBUG/%s/*/*/step*root'%(eosDir,decay))]
#inputFiles_ = ['file:%s'%path for path in glob('%s/AODSIM/%s/*/*/step*root'%(eosDir,decay))]
inputFiles_ = ['%s/%s'%(xrootd,path) for path in glob('%s/*/*root'%(eosDir))]

listname = 'list_HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4.txt'
with open(listname, 'w') as list_file:
    for inputFile in inputFiles_:
        list_file.write("%s\n" % inputFile)

maxEvents_=-1
skipEvents_=0

#cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d"%(cfg,inputFiles_,maxEvents_,skipEvents_)
cmd="cmsRun %s inputFiles_load=%s maxEvents=%d skipEvents=%d outputFile=%s/IMGs/%s_IMG.root"%(cfg,listname,maxEvents_,skipEvents_,outdir,output)
#print '%s'%cmd
os.system(cmd)

#os.system('mv cEB*.eps %s/'%(inputTag))
