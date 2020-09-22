#import pyarrow.parquet as pq
#import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
import glob, os

import dask.array as da

import h5py
from scipy.misc import imresize

import matplotlib.pyplot as plt
%matplotlib inline
from matplotlib.colors import LogNorm
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator

from skimage.measure import block_reduce
from numpy.lib.stride_tricks import as_strided

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--infile', default='output_DoubleTau.root', type=str, help='Input root file.')
args = parser.parse_args()

def loadcanvas():
  canvas = TCanvas('canvas','canvas',400,20,1400,1000)
  canvas.SetFillColor(0)
  canvas.SetBorderMode(0)
  canvas.SetFrameFillStyle(0)
  canvas.SetFrameBorderMode(0)
  canvas.SetTickx(0)
  canvas.SetTicky(0)
  return canvas

def loadlegend(top, bottom, left, right):
  relPosX    = 0.200
  relPosY    = 0.005
  posX = -1
  posX = 1 - right - relPosX*(1-left-right)
  posY = 1 - top - relPosY*(1-top-bottom)
  legendOffsetX = -0.060
  legendOffsetY = 0.
  textSize   = 0.05
  textFont   = 60
  legendSizeX = 0.5
  legendSizeY = 0.2
  legend = TLegend(posX-legendSizeX+legendOffsetX,posY-legendSizeY+legendOffsetY,posX+legendOffsetX,posY+legendOffsetY)
  legend.SetTextSize(textSize)
  legend.SetLineStyle(0)
  legend.SetBorderSize(0)
  return legend

def upsample_array(x, b0, b1):
    r, c = x.shape                                    # number of rows/columns
    rs, cs = x.strides                                # row/column strides
    x = as_strided(x, (r, b0, c, b1), (rs, 0, cs, 0)) # view as a larger 4D array

    return x.reshape(r*b0, c*b1)/(b0*b1)              # create new 2D array with same total occupancy 

def resample_EE(imgECAL, factor=2):
    # EE-
    imgEEm = imgECAL[:140-85] # EE- in the first 55 rows
    imgEEm = np.pad(imgEEm, ((1,0),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEm_dn = block_reduce(imgEEm, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEm_dn_up = upsample_array(imgEEm_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor**2
    imgECAL[:140-85] = imgEEm_dn_up[1:] ## replace the old EE- rows

    # EE+
    imgEEp = imgECAL[140+85:] # EE+ in the last 55 rows
    imgEEp = np.pad(imgEEp, ((0,1),(0,0)), 'constant', constant_values=0) # for even downsampling, zero pad 55 -> 56
    imgEEp_dn = block_reduce(imgEEp, block_size=(factor, factor), func=np.sum) # downsample by summing over [factor, factor] window
    imgEEp_dn_up = upsample_array(imgEEp_dn, factor, factor)/(factor*factor) # upsample will use same values so need to correct scale by factor*factor
    imgECAL[140+85:] = imgEEp_dn_up[:-1] # replace the old EE+ rows

    return imgECAL

def crop_jet(imgECAL, iphi, ieta, jet_shape=125):
    # NOTE: jet_shape here should correspond to the one used in RHAnalyzer
    off = jet_shape//2
    iphi = int(iphi*5 + 2) # 5 EB xtals per HB tower
    ieta = int(ieta*5 + 2) # 5 EB xtals per HB tower

    # Wrap-around on left side
    if iphi < off:
        diff = off-iphi
        img_crop = np.concatenate((imgECAL[:,ieta-off:ieta+off+1,-diff:],
                                   imgECAL[:,ieta-off:ieta+off+1,:iphi+off+1]), axis=-1)
    # Wrap-around on right side
    elif 360-iphi < off:
        diff = off - (360-iphi)
        img_crop = np.concatenate((imgECAL[:,ieta-off:ieta+off+1,iphi-off:],
                                   imgECAL[:,ieta-off:ieta+off+1,:diff+1]), axis=-1)
    # Nominal case
    else:
        img_crop = imgECAL[:,ieta-off:ieta+off+1,iphi-off:iphi+off+1]
    return img_crop



rhTreeStr = args.infile 
rhTree = ROOT.TChain("fevt/RHTree")
rhTree.Add(rhTreeStr)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> Input file:",rhTreeStr
print " >> nEvts:",nEvts
outStr = 'TauImage.root' 
print " >> Output file:",outStr

##### MAIN #####

h_Xjet = ROOT.TProfile2D("TauImage", "E(i#phi,i#eta);i#phi;i#eta", 72, 0, 72, 58,-29,29)

#  hHBHE_energy = fs->make<TProfile2D>("HBHE_energy", "E(i#phi,i#eta);i#phi;i#eta", HBHE_IPHI_NUM, HBHE_IPHI_MIN-1, HBHE_IPHI_MAX, 2*HBHE_IETA_MAX_HE,-HBHE_IETA_MAX_HE,HBHE_IETA_MAX_HE );

#
#static const int EB_IPHI_MIN = EBDetId::MIN_IPHI;//1;
#static const int EB_IPHI_MAX = EBDetId::MAX_IPHI;//360;
#static const int EB_IETA_MIN = EBDetId::MIN_IETA;//1;
#static const int EB_IETA_MAX = EBDetId::MAX_IETA;//85;
#static const int EE_MIN_IX = EEDetId::IX_MIN;//1;
#static const int EE_MIN_IY = EEDetId::IY_MIN;//1;
#static const int EE_MAX_IX = EEDetId::IX_MAX;//100;
#static const int EE_MAX_IY = EEDetId::IY_MAX;//100;
#static const int EE_NC_PER_ZSIDE = EEDetId::IX_MAX*EEDetId::IY_MAX; // 100*100
#static const int HBHE_IETA_MAX_FINE = 20;
#static const int HBHE_IETA_MAX_HB = hcaldqm::constants::IETA_MAX_HB;//16;
#static const int HBHE_IETA_MIN_HB = hcaldqm::constants::IETA_MIN_HB;//1
#static const int HBHE_IETA_MAX_HE = hcaldqm::constants::IETA_MAX_HE;//29;
#static const int HBHE_IETA_MAX_EB = hcaldqm::constants::IETA_MAX_HB + 1; // 17
#static const int HBHE_IPHI_NUM = hcaldqm::constants::IPHI_NUM;//72;
#static const int HBHE_IPHI_MIN = hcaldqm::constants::IPHI_MIN;//1;
#static const int HBHE_IPHI_MAX = hcaldqm::constants::IPHI_MAX;//72;
#static const int ECAL_IETA_MAX_EXT = 140;

# Event range to process
iEvtStart = 0
iEvtEnd   = 1
#iEvtEnd   = nEvts 
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nJets = 0
data = {} # Arrays to be written to parquet should be saved to data dict
sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    # Initialize event
    rhTree.GetEntry(iEvt)

    print " .. Processing entry",iEvt
    #if iEvt % 10000 == 0:
    #    print " .. Processing entry",iEvt

    ECAL_energy = np.array(rhTree.ECAL_energy).reshape(280,360)
    ECAL_energy = resample_EE(ECAL_energy)
    HBHE_energy = np.array(rhTree.HBHE_energy).reshape(56,72)
    HBHE_energy = upsample_array(HBHE_energy, 5, 5) # (280, 360)
    TracksAtECAL_pt = np.array(rhTree.ECAL_tracksPt).reshape(280,360)
    data['X_CMSII'] = np.stack([TracksAtECAL_pt, ECAL_energy, HBHE_energy], axis=0) # (3, 280, 360)

    # Jet attributes 
    ys  = rhTree.jetIsTau
    pts = rhTree.jetPt
    m0s = rhTree.jetM
    iphis = rhTree.jetSeed_iphi
    ietas = rhTree.jetSeed_ieta
    pdgIds = rhTree.jetPdgIds
    njets = len(ys)

    for i in range(njets):

        data['y'] = ys[i]
        data['pt'] = pts[i]
        data['m0'] = m0s[i]
        data['iphi'] = iphis[i]
        data['ieta'] = ietas[i]
        data['pdgId'] = pdgIds[i]
        data['X_jet'] = crop_jet(data['X_CMSII'], data['iphi'], data['ieta']) # (3, 125, 125)
        print "ys      = ", ys[i]
        print "pt      = ", pts[i]
        print "pdgIds  = ", pdgIds[i]
        print "iphis   = ", iphis[i]
        print "ietas   = ", ietas[i]
        for j in data['X_CMSII']:
            print j
            for k in j:
                print k
            #print "X_CMSII = ", data['X_CMSII'][0]
            #print "X_jet   = ", data['X_jet'][0]
        nJets += 1
    plt.figure(1)
    plt.subplot(221)
    plt.imshow(data['X_jet'][1,:,:,0])
    plt.title("Channel 0")  # Energy
    plt.grid(True)
    plt.subplot(222)
    plt.imshow(data['X_jet'][1,:,:,1])
    plt.title("Channel 1")  # Time
    plt.grid(True)

plt.show() 

#PLOT HERE!!
#h_Xjet = fs->make<TProfile2D>("ECAL_tracksPt", "E(i#phi,i#eta);i#phi;i#eta", EB_IPHI_MAX, EB_IPHI_MIN-1, EB_IPHI_MAX, 2*ECAL_IETA_MAX_EXT, -ECAL_IETA_MAX_EXT,   ECAL_IETA_MAX_EXT );
#hECAL_tracksPt->Fill( iphi_, ieta_signed_, trackPt_ );

print " >> nJets:",nJets
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
print "========================================================"
