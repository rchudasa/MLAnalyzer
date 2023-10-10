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

decay='DYToTauTau_M-50_13TeV'

cluster = 'FNAL'
#xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
#local='/eos/uscms/store/group/lpcml/rchudasa/NTuples/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_ntuples/230327_062100/0000'
local =''
outDir=''
if(cluster=='CERN'):
    local='/eos/cms/store/group/phys_heavyions/rchudasa/e2e/RHAnalyzer_Ntuples/DYToTauTau_M-50_13TeV-powheg_pythia8'
    outDir='/eos/user/r/rchudasa/e2e_project/ParquetFiles/DYToTauTau_M-50_13TeV-powheg_pythia8'
if(cluster=='FNAL'):
    local='/eos/uscms/store/group/lpcml/rchudasa/NTuples/DYToTauTau_M-50_13TeV-powheg_pythia8/DYToTauTau_ntuples_miniAOD/231005_120502/0000'
    outDir='/eos/uscms/store/group/lpcml/rchudasa/ParquetFiles/DYToTauTau_M-50_13TeV-powheg_pythia8/miniAODJets'


#Paths to input files
#rhFileList = '%s/%s/*.root'%(xrootd,local)
rhFileList = '%s/output*.root'%(local)
#rhFileList = '%s/*.root'%(local)
print " >> Input file list: %s"%rhFileList
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print " >> %d files found"%len(rhFileList)
#rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
#print(' >> Input File[0]: %s'%rhFileList[0])
sort_nicely(rhFileList)

#files_per_run =1 
files_per_run =10 

if not os.path.isdir(outDir):
    os.makedirs(outDir)
print ' >> Output directory: %s'%outDir

files_per_run = 10
file_idx_ = range( 0, len(rhFileList), files_per_run )
n_iter_ = len( file_idx_ )
file_idx_.append( len(rhFileList) )
print ( file_idx_ )

for irun_ in range( n_iter_ ):
    #to do run 10,42
    #if (irun_ < 8) : continue
    if (irun_ > 2) : continue
  
    print("file idx[irun]", irun_, " ", file_idx_[irun_])
    print("file idx[irun+1]", irun_+1, " ", file_idx_[irun_+1])
    files_ = rhFileList[ file_idx_[ irun_ ] : file_idx_[irun_+1] ]  
    for idx_, file_ in enumerate(files_):
        print(' >> Input File[%d]: %s' % ( idx_, file_ ) )
    
    proc_file = 'convert_root2pq_tau_jet.py'
    processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhFile, outDir, decay, ( irun_*files_per_run + i + 1 )) for i,rhFile in enumerate( files_ )]
    print(' >> Process[0]: %s'%processes[0])
    
    pool = Pool(processes=len(processes))
    pool.map(run_process, processes)

'''
proc_file = 'convert_root2pq_tau_jet.py'
#processes = []

for it,i in enumerate(range(0, len(rhFileList), files_per_run)):
    #if(it>1): continue
    #if(it>2): continue
    rhFileList_batch = rhFileList[i:i+files_per_run]
    #processes.append('%s -i %s -o %s -d %s -n %d'%(proc_file, ' '.join(rhFileList_batch), outDir, decay, it))
    processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, ' '.join(rhFileList_batch), outDir, decay, it)]
    #print(processes[it])
    #print(' >> Process[0]: %s'%processes[0])
    
    pool = Pool(processes=2)
    #pool = Pool(processes=len(processes))
    pool.map(run_process, processes)
'''
