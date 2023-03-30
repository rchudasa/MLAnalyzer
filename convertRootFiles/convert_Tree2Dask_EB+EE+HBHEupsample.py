import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da

eosDir='/eos/uscms/store/user/mba2012/IMGs'
#decays = ["H125GGgluonfusion_Pt25_Eta23_13TeV_TuneCUETP8M1_HighLumiPileUpv2", "PromptDiPhoton_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_HighLumiPileUp"]
decays = ['H125GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv4','GJet_DoubleEMEnriched_PtHat20_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUp']
#decays = ['dummy','GJet_DoubleEMEnriched_PtHat20_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_HighLumiPileUp']

#chunk_size = 250
chunk_size = 100
#scale = [100., 150.]
scale = [100., 100.]

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
def load_X_upsampled(tree, start_, stop_, branches_, readouts, scale, upscale):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_) 
    # Convert the object array X to a multidim array:
    # 1: for each event x in X, concatenate the object columns (branches) into a flat array of shape (readouts*branches)
    # 2: reshape the flat array into a stacked array: (branches, readouts)
    # 3: embed each stacked array as a single row entry in a list via list comprehension
    # 4: convert this list into an array with shape (events, branches, readouts) 
    X = np.array([np.concatenate(x).reshape(len(branches_),readouts[0]*readouts[1]) for x in X])
    #print "X.shape:",X.shape
    X = X.reshape((-1,len(branches_),readouts[0],readouts[1]))

    #print "unsampled.shape",X.shape
    X = np.stack([tile_stacked_array(x, upscale) for x in X])
    #print "upsampled.shape",X.shape
    X = np.transpose(X, [0,2,3,1])

    # Rescale
    X /= scale 
    return X

from numpy.lib.stride_tricks import as_strided

def tile_stacked_array(X, upscale):
    #print "un-tile_stacked.shape",X.shape
    X = np.stack([tile_array(x, upscale, upscale) for x in X])
    #print "tile_stacked.shape",X.shape
    return X
    
def tile_array(x, b0, b1):
    r, c = x.shape                                    # number of rows/columns
    rs, cs = x.strides                                # row/column strides 
    x = as_strided(x, (r, b0, c, b1), (rs, 0, cs, 0)) # view a as larger 4D array
    return x.reshape(r*b0, c*b1)                      # create new 2D array

for j,decay in enumerate(decays):

    if j == 0:
        pass
        #continue

    #tfile_str = '%s/%s_IMG.root'%(eosDir,decay)
    tfile_str = '%s/%s_FEVTDEBUG_IMG.root'%(eosDir,decay)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    #neff = (nevts//1000)*1000
    #neff = (nevts//100)*100
    neff = 29900 
    #neff = 100
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    # EB
    readouts = [170,360]
    branches = ["EB_energy"]
    X_EB = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[0]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EB.shape

    # EE-
    readouts = [100,100]
    branches = ["EEm_energy"]
    X_EEm = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EEm.shape

    # EE+
    readouts = [100,100]
    branches = ["EEp_energy"]
    X_EEp = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_EEp.shape

    # HBHE
    readouts = [34,72]
    branches = ["HBHE_energy_EB"]
    X_HBHE = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_HBHE.shape

    # HBHE upsample
    readouts = [34,72]
    branches = ["HBHE_energy_EB"]
    upscale = 5
    X_HBHE_up = da.concatenate([\
                da.from_delayed(\
                    load_X_upsampled(tree,i,i+chunk_size, branches, readouts, scale[1], upscale),\
                    shape=(chunk_size, readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_HBHE_up.shape

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(X_EB.shape[0], label, dtype=np.float32),\
            chunks=(chunk_size,))

    #file_out_str = "test.hdf5"
    file_out_str = "%s/%s_IMG_EBEEHBup_RH%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),neff//1000.)
    #file_out_str = "%s/%s_IMG_RH%d-%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),int(scale[1]),neff//1000.)
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X_EB': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, 'X_HBHE': X_HBHE, '/y': y}, compression='lzf')
    #da.to_hdf5(file_out_str, {'/X': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, 'X_HBHE': X_HBHE, '/y': y}, compression='lzf')
    da.to_hdf5(file_out_str, {'/X_EB': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, 'X_HBHE': X_HBHE, 'X_HBHE_up': X_HBHE_up, '/y': y}, compression='lzf')

    print " >> Done.\n"
