import numpy as np
import ROOT
from root_numpy import root2array, tree2array
#from root_pandas.readwrite import convert_to_dataframe
from dask.delayed import delayed
import dask.array as da
import dask.dataframe as df

eosDir='/eos/cms/store/user/mandrews/ML/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs_RAW'
decays = [
    #"SinglePhotonPt50",
    "SingleElectronPt50"
    #"SinglePositronPt50"
    ]
s = 32
crop_size = int(s*s)
chunk_size = 5000
n_channels = 2

@delayed
def load_X(tree, start_, stop_):
    global crop_size, s
    #branches_ = ['EB_energy','EB_time','EB_adc0','EB_adc1','EB_adc2','EB_adc3','EB_adc4','EB_adc5','EB_adc6','EB_adc7','EB_adc8','EB_adc9',]
    branches_ = ['EB_energy','EB_time']
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    X = np.array([np.concatenate(x).reshape(-1,crop_size) for x in X]) # converts the list of list to multidim array
    #X = X.reshape(-1,n_channels,crop_size)
    #X = X[:,4:,:]
    return X

input_shape = (chunk_size,n_channels,crop_size)

for j,decay in enumerate(decays):

    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_CROPS32.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_pT_CROPS32.root'%(eosDir,decay)
    tfile_str = '%s/%s_FEVTDEBUG_n125k_IMG_CROPS32.root'%(eosDir,decay)
    print " >> Opening", tfile_str
    tfile = ROOT.TFile(tfile_str)
    #tree = tfile.Get('RHTree')
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    
    neff = (nevts//chunk_size)*chunk_size
    #neff = 249000 
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    X = da.concatenate([da.from_delayed(\
            load_X(tree,i,i+chunk_size),\
            shape=input_shape,\
            dtype=np.float32) \
            for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X.shape

    label = 1
    print " >> Class label:",label
    y = da.from_array(np.full(X.shape[0], label, dtype=np.float32), chunks=(chunk_size,))

    file_out_str = "%s/%s_IMGCROPS_n%dk_RH.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk_DIGI.hdf5"%(eosDir,decay,neff//1000)
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y}, chunks=(chunk_size,s,s,2), compression='lzf')
    da.to_hdf5(file_out_str, {'/X': X, '/y': y}, compression='lzf')

    print " >> Done.\n"
