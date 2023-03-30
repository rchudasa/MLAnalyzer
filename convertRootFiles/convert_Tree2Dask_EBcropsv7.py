import numpy as np
import ROOT
from root_numpy import tree2array, root2array
from dask.delayed import delayed
from convert_Tree2Dask_utils import *
import dask.array as da
import glob

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-l', '--label', default=0, type=int, help='Decay label.')
parser.add_argument('-n', '--file_idx_start', default=1, type=int, help='File index start.')
args = parser.parse_args()

eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'
outDir='~lpcml/nobackup/mandrews' # NOTE: Space here is limited, transfer files to EOS after processing
xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
decays = [
#'DoublePi0Pt15To100_m0To1600_pythia8_noPU'
'DoublePi0Pt15To100_m0To1600_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU'
#'DoublePhotonPt50To60_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU',
#'DoublePi0Pt50To60_m000_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU',
#'DoublePi0Pt50To60_m0To1600_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU'
#'DoublePhotonPt50To60_r9gt07_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU',
#'DoublePi0Pt50To60_m000_r9gt07_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU',
#'DoublePi0Pt50To60_m0To1600_r9gt07_pythia8_2016_25ns_Moriond17MC_PoissonOOTPU'
]
#eos_date='190204_234212'
eos_date='190207_182258'

neffs = [1] 
chunk_size = 500
scale = 1.

def get_weight(m0, m0_edges, lhood):
    # m0
    if m0 >= m0_edges[-1]:
        idx_m0 = len(m0_edges)-1
    else:
        idx_m0 = np.argmax(m0 < m0_edges)-1
    return lhood[idx_m0]

def get_weight_2d(m0, pt, m0_edges, pt_edges, lhood):
    # m0
    if m0 >= m0_edges[-1]:
        idx_m0 = len(m0_edges)-1
    else:
        idx_m0 = np.argmax(m0 < m0_edges)-1
    # pt
    if pt >= pt_edges[-1]:
        idx_pt = len(pt_edges)-1
    else:
        idx_pt = np.argmax(pt < pt_edges)-1
    return lhood[idx_m0, idx_pt]

# Loop over decays
for d, decay in enumerate(decays):

    if d != args.label:
        pass
        continue
    print '>> Doing decay[%d]: %s'%(d, decay)

    #tfile_idxs = glob.glob('%s/%s*_IMG/*/*/output_*.root'%(eosDir,decay))
    #tfile_idxs = glob.glob('%s/%s*_IMG/*/*/output_1.root'%(eosDir,decay))
    tfile_idxs = glob.glob('%s/%s*_IMG/%s/0000/output_1.root'%(eosDir,decay,eos_date))
    tfile_idxs = [s.replace('.root','').split('_')[-1] for s in tfile_idxs]
    tfile_idxs = [int(i) for i in tfile_idxs]
    tfile_idxs.sort()
    #tfile_idxs = [1] # DEBUG mode: for single, local file
    print '>> File idxs:', tfile_idxs

    # Loop over root ntuples
    for n in tfile_idxs:

        if n < args.file_idx_start:
            continue

        #tfile_str = glob.glob('%s/%s*_IMG/*/*/output_%d.root'%(eosDir,decay,n))
        tfile_str = glob.glob('%s/%s*_IMG/%s/0000/output_1.root'%(eosDir,decay,eos_date))
        assert len(tfile_str) == 1, "More than 1 file of same name found in different dirs: %s"%tfile_str
        tfile_str = '%s/%s'%(xrootd,tfile_str[0])
        print " >> For input file:", tfile_str
        tfile = ROOT.TFile(tfile_str)
        tree = tfile.Get('fevt/RHTree')
        nevts = tree.GetEntries()

        #neff = (nevts//1000)*1000
        #neff = (nevts//100)*100
        #neff = 200
        neff = int(nevts)
        if neff < chunk_size:
            chunk_size = neff
        if neff > nevts:
            neff = int(nevts)
        proc_range = range(0, neff, chunk_size)
        print " >> Total events:", nevts
        print " >> Effective events:", neff

        # EB
        readouts = [170,360]
        branches = ["EB_energy"]
        X = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", X.shape

        # SC0
        readouts = [32,32]
        branches = ["SC_energy0"]
        X_crop0 = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", X_crop0.shape

        # SC1
        readouts = [32,32]
        branches = ["SC_energy1"]
        X_crop1 = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", X_crop1.shape
        X_crop0 = da.concatenate([X_crop0, X_crop1], axis=0)

        # SC0
        readouts = [32,32]
        branches = ["SC_energyT0", "SC_energyZ0"]
        X_crop_stack0 = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", X_crop_stack0.shape

        # SC1
        readouts = [32,32]
        branches = ["SC_energyT1", "SC_energyZ1"]
        X_crop_stack1 = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", X_crop_stack1.shape
        X_crop_stack0 = da.concatenate([X_crop_stack0, X_crop_stack1], axis=0)

        # SC_mass0 
        branches = ["SC_mass0"]
        y_mass0 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", y_mass0.shape

        # SC_pT0 
        branches = ["SC_pT0"]
        y_pT0 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", y_pT0.shape

        # SC_mass1 
        branches = ["SC_mass1"]
        y_mass1 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", y_mass1.shape
        y_mass0 = da.concatenate([y_mass0, y_mass1], axis=0)

        # SC_pT1 
        branches = ["SC_pT1"]
        y_pT1 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", y_pT1.shape
        y_pT0 = da.concatenate([y_pT0, y_pT1], axis=0)

        ## Likelihood weights
        #nbins = 12
        #if j == 2:
        #  h_, xs, ys = np.histogram2d(y_mass0.compute(), y_pT0.compute(), bins=nbins, range=([0.,1.6], [50., 60.])) 
        #  #print(h_)
        #  h = 1.*h_/h_.sum()
        #  #print(h)
        #  lhood = 1./h
        #  lhood = lhood/(nbins*nbins) # ensures sum_massBin_i(h*lhood) = h.sum()
        #  #print(lhood)
        #  print('sum(h_norm*lhood):',(1.*h*lhood).sum())
        #  #wgt = da.from_array(np.array([get_weight(m, xs[:-1], lhood) for m in y_mass0.compute()]), chunks=(get_chunk_size(i,neff,chunk_size),))
        #  wgt = da.from_array(np.array([get_weight_2d(m, pt, xs[:-1], ys[:-1], lhood) for m,pt in zip(y_mass0.compute(),y_pT0.compute())]), chunks=(get_chunk_size(i,neff,chunk_size),))
        #else:
        #  h_, xs = np.histogram(y_pT0.compute(), bins=nbins, range=[50., 60.]) 
        #  h = 1.*h_/h_.sum()
        #  #print(h)
        #  lhood = 1./h
        #  lhood = lhood/nbins # ensures sum_massBin_i(h*lhood) = h.sum()
        #  #print(lhood)
        #  print('sum(h_norm*lhood):',(1.*h*lhood).sum())
        #  wgt = da.from_array(np.array([get_weight(pt, xs[:-1], lhood) for pt in y_pT0.compute()]), chunks=(get_chunk_size(i,neff,chunk_size),))
        #  #wgt = da.from_array(np.ones_like(y_mass0), chunks=(get_chunk_size(i,neff,chunk_size),))*1.652721

        # SC_DR0 
        branches = ["SC_DR0"]
        y_DR0 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", y_DR0.shape

        ## SC1
        #readouts = [32,32]
        #branches = ["SC_energy1"]
        #X_crop1 = da.concatenate([\
        #            da.from_delayed(\
        #                load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale),\
        #                shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
        #                dtype=np.float32)\
        #            for i in proc_range])
        #print " >> Expected shape:", X_crop1.shape

        # pho_pT0 
        branches = ["pho_pT0"]
        pho_pT0 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", pho_pT0.shape

        # pho_E0 
        branches = ["pho_E0"]
        pho_E0 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", pho_E0.shape

        # pho_eta0 
        branches = ["pho_eta0"]
        pho_eta0 = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> Expected shape:", pho_eta0.shape

        # eventId
        branches = ["eventId"]
        eventId = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.int32)\
                    for i in proc_range])
        print " >> Expected shape:", eventId.shape

        ## Kinematics
        #branches = ["pho_pT", "pho_E", "pho_eta", "pho_phi"]
        #X_p4 = da.concatenate([\
        #            da.from_delayed(\
        #                load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
        #                shape=(get_chunk_size(i,neff,chunk_size),len(branches)),\
        #                dtype=np.float32)\
        #            for i in proc_range])
        #print " >> Expected shape:", X_p4.shape

        # Class label
        label = d
        print " >> Class label:",label
        y = da.from_array(\
                np.full(X.shape[0], label, dtype=np.float32),\
                chunks=(get_chunk_size(i,neff,chunk_size),))

        #file_out_str = "%s/%s_IMG_RH%d_n%dk_label%d.hdf5"%(eosDir,decay,int(scale),neff//1000.,label)
        #file_out_str = "%s/%s_IMGcropV4_RH%d_n%dkx2_wgt.hdf5"%(eosDir,decay,int(scale),neff//1000.)
        #file_out_str = "%s/%s_IMG/%s_IMG_RH%d_n%d_%d.hdf5"%(eosDir,decay,decay,int(scale),neff*2,n)
        file_out_str = "%s_IMG_RH%d_n%d_%d.hdf5"%(decay,int(scale),neff*2,n)
        #file_out_str = "test.hdf5"
        print " >> Writing to:", file_out_str
        #da.to_hdf5(file_out_str, {'/X': X, '/y': y, 'eventId': eventId, 'X_crop0': X_crop0, 'X_crop1': X_crop1}, compression='lzf')
        da.to_hdf5(file_out_str, {
                                  '/X': X,
                                  '/y': y,
                                  #'eventId': eventId,
                                  'X_crop0': X_crop0,
                                  'X_crop_stack0': X_crop_stack0,
                                  #'X_crop1': X_crop1
                                  #'X_p4': X_p4
                                  'y_mass': y_mass0,
                                  'y_pT': y_pT0,
                                  'y_DR': y_DR0,
                                  #'pho_pT0': pho_pT0,
                                  #'pho_E0': pho_E0,
                                  #'pho_eta0': pho_eta0
                                  #'wgt': wgt
                                  }, compression='lzf')

        print " >> Done.\n"
