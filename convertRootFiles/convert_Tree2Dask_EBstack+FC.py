import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da
from skimage.measure import block_reduce

#eosDir='/eos/uscms/store/user/mba2012/IMGs'
#eosDir='~/work/MLHEP/CMSSW_8_0_26_patch1/src/ggAnalysis/ggNtuplizer/test'
eosDir='/eos/cms/store/user/mandrews/OPENDATA/IMGs/MGG90_Eta14v2'
decays = ['DiPhotonBorn_MGG90_Eta14v2', 'GluGluHToGG_MGG90_Eta14v2', 'GJet_MGG90_Eta14v2']
#decays = ['dummy']

#chunk_size_ = 250
chunk_size_ = 200
#scale = [100., 150.]
scale = [1., 1.]

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

def block_resample_EE(X):

    return np.array([resample_EE(x) for x in X])

def resample_EE(imgECAL, factor=2):
    
    imgECAL = np.squeeze(imgECAL)
    #print('imgECAL.shape:',imgECAL.shape)
    
    # EE-
    imgEEm = imgECAL[:140-85] # EE- in the first 55 rows
    imgEEm = np.pad(imgEEm, ((1,0),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEm_dn = block_reduce(imgEEm, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEm_dn_up = tile_array(imgEEm_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor**2
    imgECAL[:140-85] = imgEEm_dn_up[1:] ## replace the old EE- rows
    
    # EE+
    imgEEp = imgECAL[140+85:] # EE+ in the last 55 rows
    imgEEp = np.pad(imgEEp, ((0,1),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEp_dn = block_reduce(imgEEp, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEp_dn_up = tile_array(imgEEp_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor*factor
    imgECAL[140+85:] = imgEEp_dn_up[:-1] # replace the old EE+ rows
    
    return np.expand_dims(imgECAL, -1)

for j,decay in enumerate(decays):

    if j != 2:
        pass
        continue

    #tfile_str = '~/work/OPENDATA/CMSSW_5_3_32/src/MLAnalyzer/output.root'
    tfile_str = '%s/%s_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/ggtree_mc_single.root'%(eosDir)
    #tfile_str = '%s/ggtree_mc.root'%(eosDir)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    #tree = tfile.Get('ggNtuplizer/EventTree')
    nevts = tree.GetEntries()
    #neff = (nevts//1000)*1000
    #neff = (nevts//100)*100
    #neff = 63000
    neff = 57400 
    #neff = int(nevts)
    chunk_size = chunk_size_
    if neff < chunk_size:
      chunk_size = neff
    if neff > nevts:
        neff = int(nevts)
        chunk_size = int(nevts)
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    # eventId
    #branches = ["event"]
    branches = ["eventId"]
    eventId = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.int32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],eventId.shape)

    # EB
    readouts = [170,360]
    branches = ["TracksPt_EB","EB_energy"]
    #branches = ["EB_energy"]
    X_EB = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale[0]),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s: %s"%(branches[0],X_EB.shape)

    # HB_EB upsample
    readouts = [34,72]
    branches = ["HBHE_energy_EB"]
    upscale = 5
    X_HBHE_EB_up = da.concatenate([\
                da.from_delayed(\
                    load_X_upsampled(tree,i,i+chunk_size, branches, readouts, scale[1], upscale),\
                    shape=(chunk_size, readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> %s(upsampled): %s"%(branches[0],X_HBHE_EB_up.shape)

    X_EB = da.concatenate([X_EB, X_HBHE_EB_up], axis=-1)
    print " >> %s: %s"%('X_EB', X_EB.shape)

    # runId
    branches = ["runId"]
    runId = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.int32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", runId.shape

    # m0
    branches = ["m0"]
    m0 = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", m0.shape

    # diPhoE
    branches = ["diPhoE"]
    diPhoE = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", diPhoE.shape

    # diPhoPt
    branches = ["diPhoPt"]
    diPhoPt = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", diPhoPt.shape


    # FC inputs
    branches = ["FC_inputs"]
    X_FC = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,5),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X_FC.shape

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(len(eventId), label, dtype=np.float32),\
            chunks=(chunk_size,))

    #file_out_str = "test%d.hdf5"%label
    file_out_str = "%s/%s_IMGall_RH%d_n%d_label%d.hdf5"%(eosDir,decay,int(scale[0]),neff,label)
    #file_out_str = "%s/%s_IMG_RH%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),neff//1000.)
    #file_out_str = "%s/%s_IMG_EBEEHBup_RH%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),neff//1000.)
    #file_out_str = "%s/%s_IMG_RH%d-%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),int(scale[1]),neff//1000.)
    print " >> Writing to:", file_out_str
    #da.to_hdf5(file_out_str, {'/X': X_EB, '/y': y, 'runId': runId, 'm0': m0, 'diPhoE': diPhoE, 'diPhoPt': diPhoPt}, compression='lzf')
    da.to_hdf5(file_out_str, {
                              '/X': X_EB,
                              '/X_FC': X_FC,
                              '/y': y,
                              'runId': runId,
                              'm0': m0,
                              'diPhoE': diPhoE,
                              'diPhoPt': diPhoPt}, compression='lzf')
    print " >> Done.\n"
