from torch.utils.data import *
from sklearn.metrics import roc_curve, auc

import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
np.random.seed(0)
import glob, os

import dask.array as da

from scipy.misc import imresize

import matplotlib.pyplot as plt
#%matplotlib inline
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec

from skimage.measure import block_reduce
from numpy.lib.stride_tricks import as_strided

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--infile', default='output_DoubleTau.root', type=str, help='Input root file.')
parser.add_argument('-n', '--nEvents', default=10, type=int, help='Number of events.')
args = parser.parse_args()

fileStr = args.infile
f0s = glob.glob(fileStr)

class ParquetDataset(Dataset):
    def __init__(self, filename):
        self.parquet = pq.ParquetFile(filename)
        self.cols = None # read all columns
        #self.cols = ['X_jet.list.item.list.item.list.item','y'] 
    def __getitem__(self, index):
        data = self.parquet.read_row_group(index, columns=self.cols).to_pydict()
        data['X_jet'] = np.float32(data['X_jet'][0])
        data['y'] = np.float32(data['y'])
        data['m0'] = np.float32(data['m0'])
        data['pt'] = np.float32(data['pt'])
        # Preprocessing
        data['X_jet'][data['X_jet'] < 1.e-3] = 0. # Zero-Suppression
        data['X_jet'][-1,...] = 25.*data['X_jet'][-1,...] # For HCAL: to match pixel intensity dist of other layers
        data['X_jet'] = data['X_jet']/100. # To standardize
        return dict(data)
    def __len__(self):
        return self.parquet.num_row_groups

def plotJet(img, mins, maxs, str_):
    plt.imshow(np.zeros_like(img[6,:,:]), cmap='Greys', vmin=0., vmax=1., alpha=0.9)
    if maxs[-1] > 0 : plt.imshow(img[6,:,:], cmap='Greens', norm=LogNorm(), alpha=0.9, vmin=mins[-1], vmax=maxs[-1])
    if maxs[-2] > 0 : plt.imshow(img[5,:,:], cmap='Greens', norm=LogNorm(), alpha=0.9, vmin=mins[-2], vmax=maxs[-2])
    if maxs[-3] > 0 : plt.imshow(img[4,:,:], cmap='Greens', norm=LogNorm(), alpha=0.9, vmin=mins[-3], vmax=maxs[-3])
    if maxs[-4] > 0 : plt.imshow(img[3,:,:], cmap='Greens', norm=LogNorm(), alpha=0.9, vmin=mins[-4], vmax=maxs[-4])
    if maxs[-5] > 0 : plt.imshow(img[2,:,:], cmap='Greys',  norm=LogNorm(), alpha=0.9, vmin=mins[-5], vmax=maxs[-5])
    if maxs[-6] > 0 : plt.imshow(img[1,:,:], cmap='Blues',  norm=LogNorm(), alpha=0.9, vmin=mins[-6], vmax=maxs[-6])
    if maxs[-7] > 0 : plt.imshow(img[0,:,:], cmap='Oranges',norm=LogNorm(), alpha=0.9, vmin=mins[-7], vmax=maxs[-7])
    #plt.colorbar(fraction=0.046, pad=0.04)
    ax = plt.axes()
    plt.xlim([0., 125.+0.])
    plt.xticks(np.arange(0,150,25))
    plt.xlabel(r"$\mathrm{i\varphi}'$", size=28) #28, 30
    ax.xaxis.set_tick_params(direction='in', which='major', length=6.)
    plt.ylim([0., 125.+0.])
    plt.yticks(np.arange(0,150,25))
    plt.ylabel(r"$\mathrm{i\eta}'$", size=28) #28, 30
    ax.yaxis.set_tick_params(direction='in', which='major', length=6.)
    #plt.savefig(str_, bbox_inches='tight')
    plt.savefig(str_, bbox_inches='tight', format='png')
    plt.show()


def plotJet_chnl(img, cmap_, xmin, xmax, str_):
#def plotJet_chnl(img, cmap_):
    plt.imshow(np.zeros_like(img), cmap='Greys', vmin=0., vmax=1., alpha=0.9)
    plt.imshow(img, cmap=cmap_, norm=LogNorm(), alpha=0.9, vmin=xmin, vmax=xmax)
    #plt.imshow(img, cmap=cmap_, norm=LogNorm(), alpha=0.9)
    ax = plt.axes()
    plt.xlim([0., 125.+0.])
    plt.xticks(np.arange(0,150,25))
    plt.xlabel(r"$\mathrm{i\varphi}'$", size=28) #28, 30
    ax.xaxis.set_tick_params(direction='in', which='major', length=6.)
    plt.ylim([0., 125.+0.])
    plt.yticks(np.arange(0,150,25))
    plt.ylabel(r"$\mathrm{i\eta}'$", size=28) #28, 30
    ax.yaxis.set_tick_params(direction='in', which='major', length=6.)
    #plt.colorbar()
    #plt.savefig(str_)
    #plt.savefig(str_, bbox_inches='tight')
    plt.savefig(str_, bbox_inches='tight', format='png')
    plt.show()

dset_train = ParquetDataset(fileStr)
train_cut = 50
idxs = np.random.permutation(len(dset_train))
train_sampler = sampler.SubsetRandomSampler(idxs[:train_cut])
#train_loader = DataLoader(dataset=dset_train, batch_size=32, num_workers=0, sampler=train_sampler, pin_memory=True)
train_loader = DataLoader(dataset=dset_train, batch_size=4, num_workers=0, shuffle=False, pin_memory=True)
print dset_train
for i, data in enumerate(train_loader):
    print("Loop ", i)
    if i > args.nEvents: break
    #print data
    #print data['X_jet']
    X_train = data['X_jet']
    print X_train.shape
    print len(X_train)
    y_train = data['y']

    plt.rcParams["font.family"] = "Helvetica"
    plt.rcParams["figure.figsize"] = (5,5)
    plt.rcParams.update({'font.size': 26})
    
    cmap = ['Oranges','Blues','Greys','Greens','Greens','Greens','Greens']

    event = 0

    #xmin, xmax = {}, {}
    #img0 = {}
    #for i in range(3):
    #    print(X_train[:,i,:,:])
    #    print(da.mean(X_train[:,i,:,:], axis=0))
    #    img0[i] = da.mean(X_train[:,i,:,:], axis=0)
    #    print(img0[i].shape)
    #    xmin[i] = da.min(img0[i]).compute()
    #    xmax[i] = da.max(img0[i]).compute()
    #    print(xmin[i], xmax[i])
    
    min_ = 0.0001

    img = X_train[event,:,:,:]
    print("JET LABEL IS  ", y_train[event])
    #Selecting only taus
    if y_train[event] == 0: continue

    for j in range(7):
    #for i in [1,0,2]:
        #print (X_train[event,i,:,:], "SHAPE : ", X_train[event,i,:,:].shape)
        img_ = img[j,:,:]
        max_ = img_.max()
        if max_ == 0: continue
        print("Channel ", j, " , Max = ", max_)
        #plotJet_chnl(img_, cmap[i], 1.e-3, img_.max(), 'images/single_event%d_jet0_chnl%d.png'%(event,i))
        plotJet_chnl(img_, cmap[j], min_, max_, 'images/tau_event%d_jet0_chnl%d.png'%(i,j))
        #plotJet_chnl(img_, cmap[i], 1.e-3, img_.max(), 'images/single_event%d_jet0_chnl%d.jpg'%(event,i))

    mins = [0.0001]*7
    maxs = [X_train[event,0,:,:].max(), X_train[event,1,:,:].max(), X_train[event,2,:,:].max(), X_train[event,3,:,:].max(), X_train[event,4,:,:].max(), X_train[event,5,:,:].max(), X_train[event,6,:,:].max()]
    print "Min = ", mins, " | Max = ", maxs
    #plotJet(X_train[event,:,:,:], mins, maxs, 'images/single_event%d_jet0.png'%(event))
    plotJet(img, mins, maxs, 'images/tau_event%d_jet0.png'%(i))
    
#print len(f0s)
#assert len(f0s) == 3, len(f0s)
#d0s = [h5py.File(f) for f in f0s]
#print(list(d0s[0].keys()))

#f1s = glob.glob('%s/QCD_Pt_80_170_00001_IMGjet_n*_event1_jet0_run?.hdf5'%(indir))
#assert len(f1s) == 3, len(f1s)
#d1s = [h5py.File(f) for f in f1s]
#print(list(d1s[0].keys()))

