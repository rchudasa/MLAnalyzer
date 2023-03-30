import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da

eosDir='/eos/uscms/store/user/mba2012/IMGs'
decays = ["H125GGgluonfusion_Pt25_Eta23_13TeV_TuneCUETP8M1_HighLumiPileUpv2", "PromptDiPhoton_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_HighLumiPileUp"]

chunk_size = 250
scale = [100., 150.]

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

for j,decay in enumerate(decays):

    if j == 1:
        continue

    tfile_str = '%s/%s_FEVTDEBUG_nXXX_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_n250k_IMG.root'%(eosDir,decay)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    neff = (nevts//1000)*1000
    #neff = 250000 
    neff = 1000
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
    readouts = [40,72]
    branches = ["HBHE_energy_EB"]
    X_HBHE = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[1]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_HBHE.shape

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(X_EB.shape[0], label, dtype=np.float32),\
            chunks=(chunk_size,))

    #file_out_str = "%s/%s_IMG_RHraw_n%dk.hdf5"%(eosDir,decay,neff//1000.)
    #file_out_str = "%s/%s_IMG_RHv1_n%dk.hdf5"%(eosDir,decay,neff//1000.)
    #file_out_str = "%s/%s_IMG_RH%d_n%dk.hdf5"%(eosDir,decay,int(rescaler),neff//1000.)
    #file_out_str = "%s/%s_IMG_RH%d-%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),int(scale[1]),neff//1000.)
    file_out_str = "%s/TEST_%s_IMG_RH%d-%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),int(scale[1]),neff//1000.)
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y}, chunks=(chunk_size,s,s,2), compression='lzf')
    #da.to_hdf5(file_out_str, {'/X_EB': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, '/y': y}, compression='lzf')
    da.to_hdf5(file_out_str, {'/X_EB': X_EB, 'X_EEm': X_EEm, 'X_EEp': X_EEp, 'X_HBHE': X_HBHE, '/y': y}, compression='lzf')

    print " >> Done.\n"
