import os, glob, re
from multiprocessing import Pool

#print("ENVIRONMENT: ", os.environ)

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

decay='HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4'

#xrootd='root://cmsxrootd.fnal.gov' # FNAL
xrootd='root://eoscms.cern.ch' # CERN
local='/eos/user/d/ddicroce/HToTauTau_m3p6To15_pT20To200_ctau0To3_eta0To1p4_RHAnalyzer_gen_unbiased'
#local='/eos/user/d/ddicroce/ML/TauClassifier/HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_biased/HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_biased_v2/201027_101416/0000/'
#local='/uscms/home/ddicroce/nobackup/TauClassifier/HToTauTau_m3p6To15_pT0To200_ctau0To3_eta0To1p4_biased_v2/201027_101416/0000'

# Paths to input files
#rhFileList = '%s/%s/*.root'%(xrootd,local)
rhFileList = '%s/*.root'%(local)
print(" >> Input file list: %s"%rhFileList)
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print(" >> %d files found"%len(rhFileList))
#rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
#print(' >> Input File[0]: %s'%rhFileList[0])
sort_nicely(rhFileList)

files_per_run = 10
file_idx_ = range( 0, len(rhFileList), files_per_run )
n_iter_ = len( file_idx_ )
file_idx_.append( len(rhFileList) )
print ( file_idx_ )

for irun_ in range( n_iter_ ):
    #if (irun_ < 8) : continue
    files_ = rhFileList[ file_idx_[ irun_ ] : file_idx_[irun_+1] ]  
    for idx_, file_ in enumerate(files_):
        print(' >> Input File[%d]: %s' % ( idx_, file_ ) )
    
    # Output path
    outDir='/afs/cern.ch/user/d/ddicroce/work/ML/TauIdentifier/CMSSW_10_2_20_UL/src/MLAnalyzer/IMG'
    if not os.path.isdir(outDir):
        os.makedirs(outDir)
    print(' >> Output directory: %s'%outDir)
    
    proc_file = 'convert_root2pq_jet.py'
    processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhFile, outDir, decay, ( irun_*files_per_run + i + 1 )) for i,rhFile in enumerate( files_ )]
    print(' >> Process[0]: %s'%processes[0])
    
    pool = Pool(processes=len(processes))
    pool.map(run_process, processes)
