import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da

eosDir='/eos/uscms/store/user/mba2012/FCs/h24gamma_eta14'
#decays = ['h22gammaSM_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
#decays = ['SM2gamma_1j_1M_noPU', 'h24gamma_1j_1M_1GeV_noPU']
decays = ['SM2gamma_1j_1M_noPU', 'h22gammaSM_1j_1M_noPU']

chunk_size = 1000
scale = 1.

@delayed
def load_X(tree, start_, stop_, branches_, readouts, scale):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    # Convert the object array X to a multidim array:
    # 1: for each event x in X, concatenate the object columns (branches) into a flat array of shape (readouts*branches)
    # 2: reshape the flat array into a stacked array: (branches, readouts)
    # 3: embed each stacked array as a single row entry in a list via list comprehension
    # 4: convert this list into an array with shape (events, branches, readouts) 
    X = np.array([np.concatenate(x).reshape(len(branches_),readouts[0]*readouts[1]) for x in X])
    #print "X.shape:",X.shape
    X = X.reshape((-1,len(branches_),readouts[0],readouts[1]))
    X = np.transpose(X, [0,2,3,1])

    # Rescale
    X /= scale 
    return X

@delayed
def load_single(tree, start_, stop_, branches_):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    X = np.array([x[0] for x in X])

    return X

for j,decay in enumerate(decays):

    if j == 0:
        pass
        continue

    tfile_str = '%s/%s_FEVTDEBUG_IMG.root'%(eosDir,decay)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    neff = (nevts//1000)*1000
    #neff = 40 
    #neff = 170000
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    # FC inputs 
    branches = ["FC_inputs"]
    X = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,5),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X.shape

    # eventId
    branches = ["eventId"]
    eventId = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.int32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", eventId.shape

    # m0
    branches = ["m0"]
    m0 = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", m0.shape

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(X.shape[0], label, dtype=np.float32),\
            chunks=(chunk_size,))

    file_out_str = "%s/%s_FC_n%dk_label%d.hdf5"%(eosDir,decay,neff//1000.,label)
    #file_out_str = "test.hdf5"
    print " >> Writing to:", file_out_str
    da.to_hdf5(file_out_str, {'/X': X, '/y': y, 'eventId': eventId, 'm0': m0}, compression='lzf')

    print " >> Done.\n"
