import numpy as np
import ROOT
from root_numpy import root2array, tree2array
#from root_pandas.readwrite import convert_to_dataframe
from dask.delayed import delayed
import dask.array as da
import dask.dataframe as df
import h5py

eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs'
#eosDir='/eos/cms/store/user/mandrews/ML/IMGs_RAW'
decays = ["SinglePhotonPt50","SingleElectronPt50"]
#decays = ["DoublePhotonFlatPt10To60","DoubleElectronFlatPt10To60"]
s = 32
#crop_size = int(s*s)
crop_size = 170*360
chunk_size = 100
n_channels = 2
ver = 1

chunk_shape = (chunk_size,2,crop_size)

@delayed
def process_chunk(x):

    # Energy
    E = x[:,0,:]
    if ver == 1:
        E /= 50.
        #E /= 120.
    if ver == 3:
        nonzero = (E > 0.)
        E[nonzero] = (np.log10(E[nonzero])+4.) / 4.
    if ver == 4:
        nonzero = (E > 0.)
        E[nonzero] = (np.log10(E[nonzero])+1.) / 4.
    if ver == 5:
        norm = np.expand_dims(np.linalg.norm(E, axis=1), -1)
        E /= norm 
        E[np.isnan(E)] = 0.

    # Time
    t = x[:,1,:]
    t /= 50.

    X = np.stack([E, t], axis=1)
    X = X.reshape((-1,X.shape[1],170,360)) 
    X = np.transpose(X, [0,2,3,1]) 
    return X

for j,decay in enumerate(decays):

    #file_in_str = "%s/%s_IMGCROPS_n249k_RH.hdf5"%(eosDir,decay)
    #file_in_str = "%s/%s_IMGCROPS_n249k_pT.hdf5"%(eosDir,decay)
    file_in_str = "%s/%s_n250k_IMG_RHraw.hdf5"%(eosDir,decay)
    dset = h5py.File(file_in_str)
    X_in = da.from_array(dset['/X'], chunks=chunk_shape)
    y_in = da.from_array(dset['/y'], chunks=(chunk_size,))
    assert X_in.shape[0] == y_in.shape[0]
    events = X_in.shape[0]
    #events = 10000
    assert events % chunk_size == 0
    print " >> Doing decay:", decay
    print " >> Input file:", file_in_str
    print " >> Total events:", events

    print " >> Processing..."
    X = da.concatenate([da.from_delayed(\
                        process_chunk(X_in[i:i+chunk_size]),\
                        shape=(chunk_size,170,360,2),\
                        dtype=np.float32)\
                        for i in range(0,events,chunk_size)])

    #file_out_str = "%s/%s_IMGCROPS_n249k_pT_RHv%d.hdf5"%(eosDir,decay,ver)
    file_out_str = "%s/%s_n250k_IMG_RHv%d.hdf5"%(eosDir,decay,ver)
    print " >> Writing to:", file_out_str
    da.to_hdf5(file_out_str, {'/X': X, '/y': y_in}, compression='lzf')

    print " >> Done.\n"
