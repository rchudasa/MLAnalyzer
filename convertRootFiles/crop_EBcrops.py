import numpy as np
from os.path import splitext

import matplotlib.pyplot as plt
import ROOT
from scipy.ndimage import maximum_position
from scipy.sparse import csr_matrix
#import root_numpy
from matplotlib.colors import LogNorm
import argparse

# Register command line options
parser = argparse.ArgumentParser(description='Run STEALTH selection.')
parser.add_argument('-d','--decay', required=True, help='Decay:Single*Pt50_FEVTDEBUG_n250k_IMG',type=str)
parser.add_argument('-n','--nevts', default=10, help='Number of events to process.',type=int)
args = parser.parse_args()

def crop_around_max(b,row0_,col0_):
    global n_rows, n_cols, w
    return np.array(b, dtype=np.float32).reshape(n_rows,n_cols)[row0_-w:row0_+w,col0_-w:col0_+w].flatten()

s = 32
crop_size = int(s*s)
w = s//2
n_rows = 170 # n_phi
n_cols = 360 # n_eta

##### I/O #####
eos_dir = '/eos/cms/store/user/mandrews/ML'
#decay = 'SinglePhotonPt50_FEVTDEBUG_n250k_IMG'
#decay = 'SingleElectronPt50_FEVTDEBUG_n250k_IMG'
decay = args.decay

file_in_str = '%s/IMGs/%s.root'%(eos_dir,decay)
tree_in = ROOT.TChain('fevt/RHTree')
tree_in.Add(file_in_str)
n_events = tree_in.GetEntries()
branch_list = [br.GetName() for br in tree_in.GetListOfBranches()]
print " >> Read input file:", file_in_str
print " >> N of events:", n_events
print " >> Input branch list:",branch_list

#file_out_str = 'test.root'
file_out_str = '%s/IMGs_RAW/%s_CROPS32.root'%(eos_dir,decay)
file_out = ROOT.TFile(file_out_str, 'RECREATE')
RHTree = ROOT.TTree("RHTree", "RecHit tree")
EB_energy    = np.zeros(crop_size, dtype=np.float32)
EB_time      = np.zeros(crop_size, dtype=np.float32)
EB_adc0     = np.zeros(crop_size, dtype=np.float32)
EB_adc1     = np.zeros(crop_size, dtype=np.float32)
EB_adc2     = np.zeros(crop_size, dtype=np.float32)
EB_adc3     = np.zeros(crop_size, dtype=np.float32)
EB_adc4     = np.zeros(crop_size, dtype=np.float32)
EB_adc5     = np.zeros(crop_size, dtype=np.float32)
EB_adc6     = np.zeros(crop_size, dtype=np.float32)
EB_adc7     = np.zeros(crop_size, dtype=np.float32)
EB_adc8     = np.zeros(crop_size, dtype=np.float32)
EB_adc9     = np.zeros(crop_size, dtype=np.float32)
RHTree.Branch('EB_energy'  ,EB_energy   , 'EB_energy[%d]/F'%crop_size   )
RHTree.Branch('EB_time'    ,EB_time     , 'EB_time[%d]/F'%crop_size     )
RHTree.Branch('EB_adc0'    ,EB_adc0     , 'EB_adc0[%d]/F'%crop_size    )
RHTree.Branch('EB_adc1'    ,EB_adc1     , 'EB_adc1[%d]/F'%crop_size    )
RHTree.Branch('EB_adc2'    ,EB_adc2     , 'EB_adc2[%d]/F'%crop_size    )
RHTree.Branch('EB_adc3'    ,EB_adc3     , 'EB_adc3[%d]/F'%crop_size    )
RHTree.Branch('EB_adc4'    ,EB_adc4     , 'EB_adc4[%d]/F'%crop_size    )
RHTree.Branch('EB_adc5'    ,EB_adc5     , 'EB_adc5[%d]/F'%crop_size    )
RHTree.Branch('EB_adc6'    ,EB_adc6     , 'EB_adc6[%d]/F'%crop_size    )
RHTree.Branch('EB_adc7'    ,EB_adc7     , 'EB_adc7[%d]/F'%crop_size    )
RHTree.Branch('EB_adc8'    ,EB_adc8     , 'EB_adc8[%d]/F'%crop_size    )
RHTree.Branch('EB_adc9'    ,EB_adc9     , 'EB_adc9[%d]/F'%crop_size    )
branch_list = [br.GetName() for br in RHTree.GetListOfBranches()]
print " >> Output file:",file_out_str
print " >> Output branch list:",branch_list

# Temp arrays to calculate presample
EB_adc0_ = np.zeros(crop_size, dtype=np.float32)
EB_adc1_ = np.zeros(crop_size, dtype=np.float32)
EB_adc2_ = np.zeros(crop_size, dtype=np.float32)

##### IMAGE SELECTION #####
istart, istop = 0, n_events
#if istop < args.nevts:
#    istop = n_evts
#else:
#    istop = args.nevts
row0, col0 = -1, -1
print " >> Processing entries: [",istart,"->",istop,")"
for ievt in range(istart,istop):

    # Initialize event
    if ievt > istop:
        break
    treeStatus = tree_in.LoadTree(ievt)
    if treeStatus < 0:
        break
    evtStatus = tree_in.GetEntry(ievt)
    if evtStatus <= 0:
        continue
    if ievt % 1000 == 0:
        print " .. Processing entry",ievt
    
    ### Crop around shower max ###
    
    # Get position of max adc
    row0, col0 = -1, -1
    #print len(maximum_position(np.array(tree_in.EB_adc6, dtype=np.float32).reshape(n_rows,n_cols)))
    row0, col0 = maximum_position(np.array(tree_in.EB_adc6, dtype=np.float32).reshape(n_rows,n_cols))
    if col0 < w or col0 >= n_cols-w or row0 < w or row0 >= n_rows-w:
        continue
    #print row0,col0

    ### Energy ###
    b = crop_around_max(tree_in.EB_energy,row0,col0)
    for i,val in enumerate(b):
        EB_energy[i] = val

    ### Timing ###
    b = crop_around_max(tree_in.EB_time,row0,col0)
    for i,val in enumerate(b):
        EB_time[i] = val

    ### Digis ###
    '''
    b = crop_around_max(tree_in.EB_adc0,row0,col0)
    for i,val in enumerate(b):
        EB_adc0[i] = val
    
    b = crop_around_max(tree_in.EB_adc1,row0,col0)
    for i,val in enumerate(b):
        EB_adc1[i] = val
    
    b = crop_around_max(tree_in.EB_adc2,row0,col0)
    for i,val in enumerate(b):
        EB_adc2[i] = val

    EB_adc0_ = crop_around_max(tree_in.EB_adc0,row0,col0)
    EB_adc1_ = crop_around_max(tree_in.EB_adc1,row0,col0)
    EB_adc2_ = crop_around_max(tree_in.EB_adc2,row0,col0)

    presample = np.mean([EB_adc0_, EB_adc1_, EB_adc2_], axis=0)
    presample = log_noise(presample)
    presample = 0. 
    '''

    b = crop_around_max(tree_in.EB_adc0,row0,col0)
    for i,val in enumerate(b):
        EB_adc0[i] = val
    
    b = crop_around_max(tree_in.EB_adc1,row0,col0)
    for i,val in enumerate(b):
        EB_adc1[i] = val
    
    b = crop_around_max(tree_in.EB_adc2,row0,col0)
    for i,val in enumerate(b):
        EB_adc2[i] = val
    
    b = crop_around_max(tree_in.EB_adc3,row0,col0)
    for i,val in enumerate(b):
        EB_adc3[i] = val
    
    b = crop_around_max(tree_in.EB_adc4,row0,col0)
    for i,val in enumerate(b):
        EB_adc4[i] = val
    
    b = crop_around_max(tree_in.EB_adc5,row0,col0)
    for i,val in enumerate(b):
        EB_adc5[i] = val
    
    b = crop_around_max(tree_in.EB_adc6,row0,col0)
    for i,val in enumerate(b):
        EB_adc6[i] = val
    
    b = crop_around_max(tree_in.EB_adc7,row0,col0)
    for i,val in enumerate(b):
        EB_adc7[i] = val
    
    b = crop_around_max(tree_in.EB_adc8,row0,col0)
    for i,val in enumerate(b):
        EB_adc8[i] = val

    b = crop_around_max(tree_in.EB_adc9,row0,col0)
    for i,val in enumerate(b):
        EB_adc9[i] = val

    RHTree.Fill()
    
    '''
    # Check plots
    img = b['EB_adc9'].reshape(32,32)
    plt.imshow(img, interpolation="None", cmap='seismic', vmin=-0.1, vmax=0.8) # Blues, seismic
    plt.colorbar()
    plt.show()
    hist = ROOT.TH1F("h","h",100,-1.*(img.ravel().max()+2.),img.ravel().max()+2.)
    print img.ravel().max()
    for i in img.ravel():
        hist.Fill(i)
    c = ROOT.TCanvas("c")
    ROOT.gPad.SetLogy()
    hist.Draw()
    c.Draw()
    '''

file_out.Write()
file_out.Close()
