import os, glob, re
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

xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
local='/afs/cern.ch/user/d/ddicroce/work/ML/TauIdentifier/CMSSW_9_4_17/src/MLAnalyzer/input'

decay = 'input'

# Paths to input files
rhFileList = '%s/%s/*.root'%(local, decay)
print(" >> Input file list: %s"%rhFileList)
rhFileList = glob.glob(rhFileList)
assert len(rhFileList) > 0
print(" >> %d files found"%len(rhFileList))
#rhFileList = [('%s/%s'%(xrootd, rhFile)).replace('/eos/uscms','') for rhFile in rhFileList]
print(' >> Input File[0]: %s'%rhFileList[0])
sort_nicely(rhFileList)

# Output path
outDir='/afs/cern.ch/user/d/ddicroce/work/ML/TauIdentifier/CMSSW_9_4_17/src/MLAnalyzer/'
outDir='%s/%s'%(outDir, decay)
if not os.path.isdir(outDir):
    os.makedirs(outDir)
print(' >> Output directory: %s'%outDir)

proc_file = 'convert_root2pq_jet.py'
processes = ['%s -i %s -o %s -d %s -n %d'%(proc_file, rhFile, outDir, decay, i+1) for i,rhFile in enumerate(rhFileList)]
print(' >> Process[0]: %s'%processes[0])

pool = Pool(processes=len(processes))
pool.map(run_process, processes)
