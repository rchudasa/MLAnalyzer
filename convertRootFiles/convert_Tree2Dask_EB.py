import numpy as np
import ROOT
from root_numpy import root2array, tree2array
#from root_pandas.readwrite import convert_to_dataframe
from dask.delayed import delayed
import dask.array as da
import dask.dataframe as df

eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs_RAW'
decays = ["SinglePhotonPt50","SingleElectronPt50"]
#decays = ["DoublePhotonFlatPt10To60","DoubleElectronFlatPt10To60"]
s = 32
#crop_size = int(s*s)
crop_size = 170*360
chunk_size = 100
#n_channels = 14
n_channels = 2

'''
@delayed
def load_X(tree, start_, stop_):
    global crop_size, s
    X = tree2array(tree, start=start_, stop=stop_) 
    X = np.array([np.concatenate(x).reshape(-1,crop_size) for x in X]) # converts the list of list to multidim array
    X = X.reshape((-1,X.shape[1],s,s))
    X = np.transpose(X,(0,2,3,1))
    #X = np.swapaxes(X,1,2)
    #X = np.swapaxes(X,2,3)
    return X
'''
@delayed
def load_X(tree, start_, stop_):
    global crop_size, s
    X = tree2array(tree, start=start_, stop=stop_, branches=['EBenergyRed','EBtimeRed']) 
    #print X.dtype.names
    X = np.array([np.concatenate(x).reshape(-1,crop_size) for x in X]) # converts the list of list to multidim array
    #X = X[:,:2,:]
    return X

input_shape = (chunk_size,n_channels,crop_size)

for j,decay in enumerate(decays):

    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_CROPS32.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG_pT_CROPS32.root'%(eosDir,decay)
    tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG.root'%(eosDir,decay)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    neff = (nevts//chunk_size)*chunk_size
    #neff = 100 
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    X = da.concatenate([da.from_delayed(load_X(tree,i,i+chunk_size),shape=input_shape, dtype=np.float32) \
                        for i in range(0,neff,chunk_size)])
    #data = [df.from_pandas(get_Xy(tree,i,i+chunk_size)) for i in range(0,nevts,chunk_size)]
    #data = da.concatenate([da.from_delayed(get_Xy(tree,i,i+chunk_size, nevts), shape=(chunk_size,s,s,14), dtype=np.float32) for i in range(0,neff,chunk_size)])
    #data = da.concatenate([da.from_array(get_Xy(tree,i,i+chunk_size, nevts), (200,32,32,2)) for i in range(0,neff,chunk_size)])
    #data = da.concatenate([da.from_delayed(get_Xy(tree,i,i+chunk_size), dtype=np.float32) for i in range(0,nevts,chunk_size)])

    print " >> Expected shape:", X.shape
    print " >> Class label:",j
    y = da.from_array(np.full(X.shape[0], j, dtype=np.float32), chunks=(chunk_size,))

    file_out_str = "%s/%s_n%dk_IMG_RHraw.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk_pT.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk.hdf5"%(eosDir,decay,neff//1000)
    #file_out_str = "%s/%s_IMGCROPS_n%dk_DIGI.hdf5"%(eosDir,decay,neff//1000)
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y}, chunks=(chunk_size,s,s,2), compression='lzf')
    da.to_hdf5(file_out_str, {'/X': X, '/y': y}, compression='lzf')

    print " >> Done.\n"
