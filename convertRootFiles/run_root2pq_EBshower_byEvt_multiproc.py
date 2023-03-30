import os, glob, re
import shutil
from multiprocessing import Pool

def alphanum_key(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [int(c) if c.isdigit() else c for c in re.split('([0-9]+)',s)]

def sort_nicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphanum_key)

def run_process(process):
    os.system('python %s'%process)

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-d', '--genDR', default=10, type=int, help='gen-level dR.')
parser.add_argument('-p', '--p_drop', default=1.00, type=float, help='p(drop) scale.')
#parser.add_argument('-w', '--wgt_file', default=None, type=str, help='Weight file.')
parser.add_argument('-w', '--wgt_files', default=None, nargs='+', type=str, help='Weight file.')
args = parser.parse_args()

genDR = args.genDR
p_drop_scale = args.p_drop
wgt_files = args.wgt_files

#xrootd='root://cmsxrootd.fnal.gov' # FNAL
xrootd='root://cmseos.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'
#eosDir='/eos/uscms/store/user/lpcml/mandrews'

decay = 'h24gamma_1j_1M_100MeV_PU2017_MINIAODSIM_ext3_bdt'
decay = 'h24gamma_1j_1M_200MeV_PU2017_MINIAODSIM_ext3_bdt'
decay = 'h24gamma_1j_1M_400MeV_PU2017_MINIAODSIM_ext3_bdt'
decay = 'h24gamma_1j_1M_1GeV_PU2017_MINIAODSIM_ext3_bdt'
date_str = '190807_183309'
date_str = '190807_184410'
date_str = '190807_184430'
date_str = '190807_184451'
#decay = 'GluGluHToGG_M125_13TeV_amcatnloFXFX_pythia8_MINIAODSIM_noExt3_2presel_wrapfix'
#decay = 'DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8_MINIAODSIM_wgt'
#decay = 'DYToEE_M-50_NNPDF31_13TeV-powheg-pythia8_MINIAODSIM_noExt3_2presel_wrapfix'
#decay = 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_MINIAODSIM_wgt'
#decay = 'GJet_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_MINIAODSIM_wgt'
#decay = 'QCD_Pt-30to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_MINIAODSIM_wgt'
#decay = 'QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_MINIAODSIM_wgt'
#date_str = ''
#date_str = '190812_120427'
#date_str = ''
#date_str = '190809_200713'
#date_str = '190809_200635'
#date_str = '190809_200829'
#date_str = '190809_200741'
#decay = 'Run2017B_MINIAOD_noPtOmGG_noPresel'
#decay = 'Run2017C_MINIAOD_noPtOmGG_noPresel'
#decay = 'Run2017D_MINIAOD_noPtOmGG_noPresel'
#decay = 'Run2017E_MINIAOD_noPtOmGG_noPresel'
#decay = 'Run2017F_MINIAOD_noPtOmGG_noPresel'
#date_str = '190805_020913'
#date_str = '190805_020845'
#date_str = '190805_020747'
#date_str = '190805_020716'
#date_str = '190805_020614'
#decay = 'QCD_Pt-40toInf_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_MINIAODSIM_MCstudy'
#decay = 'GJet_Pt-20to40_DoubleEMEnriched_MGG-80toInf_TuneCP5_13TeV_Pythia8_MINIAODSIM_MCstudy'
#date_str = '190627_223251'
#date_str = '190619_045146'

# Paths to input files
rhFileList = '%s/%s/%s/*/output_*.root'%(eosDir, decay, date_str)
#rhFileList = '%s/%s/%s/%s/*/output_*.root'%(eosDir, decay.partition('_MINIAODSIM')[0], decay, date_str)
#rhFileList = '%s/DoubleEG/%s/%s/*/output_*.root'%(eosDir, decay, date_str)
#print(" >> Input file list: %s"%rhFileList)
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print(" >> %d files found"%len(rhFileList))
#rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList][:10]
#rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList][:500]
rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
#rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList][:1]
print(' >> Input File[0]: %s'%rhFileList[0])
sort_nicely(rhFileList)

# Weights file
if wgt_files is not None:
  #wgt_file = '%s_mvpt_weights_pdrop%.2f.npz'%(decay, p_drop_scale)
  #wgt_files = args.wgt_files
  for wgt_file in wgt_files:
      assert os.path.isfile(wgt_file)

# Output path
outDir='/uscms/physics_grp/lpcml/nobackup/mandrews' # NOTE: Space here is limited, transfer files to EOS after processing
outDir='%s/%s'%(outDir, decay)
if not os.path.isdir(outDir):
    os.makedirs(outDir)
print(' >> Output directory: %s'%outDir)

proc_file = 'convert_root2pq_EBshower_byEvt.py'
if wgt_files is not None:
    #processes = ['%s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFile, outDir, decay, i+1, wgt_files) for i,rhFile in enumerate(rhFileList)]
    processes = ['%s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFile, outDir, decay, i+1, ' '.join(wgt_files)) for i,rhFile in enumerate(rhFileList)]
else:
  #processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhFile, outDir, decay, i+1) for i,rhFile in enumerate(rhFileList)]
  #processes = ['%s -i %s -o %s -d %s'%(proc_file, ' '.join(rhFileList), outDir, decay)]
  processes = []
  for it,i in enumerate(range(0, len(rhFileList), 500)):
    rhFileList_batch = rhFileList[i:i+500]
    #print(it, i, i+500, len(rhFileList_batch))
    processes.append('%s -i %s -o %s -d %s -n %d'%(proc_file, ' '.join(rhFileList_batch), outDir, decay, it))
#print(' >> Process[0]: %s'%processes[0])

#os.system('python %s -i %s -o %s -d %s -n %d -w %s'%(proc_file, rhFileList[0], outDir, decay, 1, wgt_file))
pool = Pool(processes=len(processes))
pool.map(run_process, processes)
pool.close()
pool.join()

for f in glob.glob('%s/%s*.parquet.*'%(outDir, decay)):
  pass
  shutil.move(f, '%s/%s/%s/'%(eosDir, decay, date_str))
  #shutil.move(f, '%s/%s/%s/%s/'%(eosDir, decay.partition('_MINIAODSIM')[0], decay, date_str))
  #shutil.move(f, '%s/DoubleEG/%s/%s/'%(eosDir, decay, date_str))
