import numpy as np
#import ROOT
#from scipy.sparse import csr_matrix
from root_numpy import root2array, tree2array
#from root_pandas import read_root
#import tables

eosDir='/eos/cms/store/user/mandrews/ML/IMGs/case1'

decays = ['SinglePhotonPt50_FEVTDEBUG','SingleElectronPt50_FEVTDEBUG']
#decay = "SingleElectronPt50_FEVTDEBUG"
#decay = "SinglePhotonPt50_FEVTDEBUG"

def x_arr(decay,start_,stop_):
  #arr = root2array('%s/%s_n250k_IMG.root'%(eosDir,decay), 'fevt/RHTree', start=start_, stop=stop_) 
  #arr = root2array('%s/%s_n250k_IMG_CROPS32.root'%(eosDir,decay), 'RHTree', start=start_, stop=stop_) 
  arr = root2array('output_numEvent1.root','fevt/RHTree', start=start_, stop=stop_) 
  print arr.dtype.names
  return np.array([np.concatenate(x).reshape(-1,61200) for x in arr])
  #return np.array([np.concatenate(x).reshape(-1,1024) for x in arr])

n_evts = 250000
chunk_size = 25000
n_chunks = int(np.ceil(float(n_evts)/float(chunk_size)))
for d in decays:
    #for i in range(n_chunks):
    for i in range(1):
      print " >> Doing chunk",i
      x = x_arr(d,i*chunk_size,(i+1)*chunk_size)
      print x.shape
      np.savez_compressed('Data_EB.npz', x=x)
      #np.savez_compressed('IMGs/%s_n250k_IMG_%d'%(d,i),x=x)
      #np.savez_compressed('IMGs/%s_n250k_IMG_CROPS32_%d'%(d,i),x=x)
