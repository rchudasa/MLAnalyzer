import numpy as np
import dask.array as da
from dask.delayed import delayed
from root_numpy import tree2array, root2array
from skimage.measure import block_reduce

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
def load_vector(tree, start_, stop_, branches_, idx_=None):
    X = tree2array(tree, start=start_, stop=stop_, branches=branches_)
    X = np.array([np.concatenate(x) for x in X])
    if idx_ is not None:
        X = X[:,idx_]
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
#    return x.reshape(r*b0, c*b1)/(b0*b1)              # create new 2D array

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

def crop_jet(imgECAL, iphi, ieta, jet_shape):
    iphi = int(iphi*5 + 2)
    ieta = int(ieta*5 + 2)
    off = jet_shape//2
    #print('iphi:%d, ieta:%d'%(iphi, ieta))
    if iphi < off:
        #print('1')
        diff = off-iphi
        #diff = diff+150
        #print(diff)
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,-diff:],
                                   imgECAL[ieta-off:ieta+off+1,:iphi+off+1]), axis=1)
    elif 360-iphi < off:
        #print('2')
        diff = off - (360-iphi)
        #diff = diff+150
        #print('diff:',diff)
        img_crop = np.concatenate((imgECAL[ieta-off:ieta+off+1,iphi-off:], 
                                   imgECAL[ieta-off:ieta+off+1,:diff+1]), axis=1)
    else:
        #print('0')
        img_crop = imgECAL[ieta-off:ieta+off+1,iphi-off:iphi+off+1]
    #print(img_crop.shape)
    return img_crop

@delayed
def crop_jet_block(Xs, iphis, ietas, jet_shape):
    X = np.array([crop_jet(x,iphi,ieta,jet_shape) for x,iphi,ieta in zip(Xs,iphis,ietas)])
    #print X.shape
    return X

def get_chunk_size(i, neff, chunk_size):
    # If last chunk less than chunk_size
    if neff-i < chunk_size:
        eff_chunk_size = neff-i
    # Otherwise, use chunk_size
    else:
        eff_chunk_size = chunk_size
    return eff_chunk_size
