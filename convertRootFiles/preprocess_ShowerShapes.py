import numpy as np
from scipy.ndimage import maximum_position, center_of_mass
from dask.delayed import delayed
import dask.array as da
import dask.dataframe as df
import h5py

#eosDir='/eos/cms/store/user/mandrews/ML/IMGs'
eosDir='/eos/cms/store/user/mandrews/ML/IMGs_RAW'
decays = [
    "SinglePhotonPt50",
    "SingleElectronPt50"
    ]
s = 32
crop_size = int(s*s)
chunk_size = 9960
#chunk_size = 19920
#chunk_size = 1000
n_channels = 4

chunk_shape = (chunk_size,n_channels,crop_size)

@delayed
def process_chunk(x):

    # Energy
    E = x[:,2,:]
    E = E.reshape((-1,s,s))
    #print E.shape
    e3x3 = np.zeros(E.shape[0], dtype=np.float32)
    e5x5 = np.zeros(E.shape[0], dtype=np.float32)
    r2x5 = np.zeros(E.shape[0], dtype=np.float32)
    r9 = np.zeros(E.shape[0], dtype=np.float32)
    etaWidth = np.zeros(E.shape[0], dtype=np.float32)
    phiWidth = np.zeros(E.shape[0], dtype=np.float32)
    eSC = np.zeros(E.shape[0], dtype=np.float32)
    eMax = np.zeros(E.shape[0], dtype=np.float32)
    for i,e in enumerate(E):
        r, c = maximum_position(e)
        e2x5_ = e[r-2:r+3,c:c+2].sum()
        e3x3_ = e[r-1:r+2,c-1:c+2].sum()
        e5x5_ = e[r-2:r+3,c-2:c+3].sum()

        eSC_ = e.sum()
        eMax_ = e[r,c].flatten()
        etaSC, phiSC = center_of_mass(e)
        #print etaSC, phiSC, eSC_
        etaWidth_, phiWidth_ = 0., 0.
        for idx, en in np.ndenumerate(e):
            #print idx[0], idx[1], en
            etaWidth_ += en/eSC_ * ((idx[0]-etaSC)*0.0174)**2
            phiWidth_ += en/eSC_ * ((idx[1]-phiSC)*0.0174)**2

        # Write
        e3x3[i] = e3x3_
        e5x5[i] = e5x5_
        r2x5[i] = e2x5_/e5x5_
        r9[i] = e3x3_/eSC_
        etaWidth[i] = np.sqrt(etaWidth_)
        phiWidth[i] = np.sqrt(phiWidth_)
        eSC[i] = eSC_
        eMax[i] = eMax_

    X = np.stack([e3x3, e5x5, r2x5, r9, etaWidth, phiWidth, eSC, eMax], axis=1)
    return X

for j,decay in enumerate(decays):

    file_in_str = "%s/%s_IMGCROPS_n249k_RH.hdf5"%(eosDir,decay)
    dset = h5py.File(file_in_str)
    X_in = da.from_array(dset['/X'], chunks=chunk_shape)
    y_in = da.from_array(dset['/y'], chunks=(chunk_size,))
    #X_in = dset['/X'][:chunk_size,...]
    #y_in = dset['/y'][:chunk_size]
    assert X_in.shape[0] == y_in.shape[0]
    events = X_in.shape[0]
    #events = 1000
    assert events % chunk_size == 0
    print " >> Doing decay:", decay
    print " >> Input file:", file_in_str
    print " >> Total events:", events

    print " >> Processing..."
    #X = np.concatenate([process_chunk(X_in[i:i+chunk_size]) for i in range(0,events,chunk_size)])
    X = da.concatenate([da.from_delayed(\
                        process_chunk(X_in[i:i+chunk_size]),\
                        shape=(chunk_size,8),\
                        dtype=np.float32)\
                        for i in range(0,events,chunk_size)])

    file_out_str = "%s/%s_PID_n249k.hdf5"%(eosDir,decay)
    print " >> Writing to:", file_out_str
    da.to_hdf5(file_out_str, {'/X': X, '/y': y_in}, compression='lzf')

    print " >> Done.\n"
