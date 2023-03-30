import glob, os
import ROOT
import numpy as np
import dask.array as da
from convert_Tree2Dask_utils import *

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-l', '--label', required=True, type=int, help='Decay label.')
parser.add_argument('-n', '--file_idx_start', default=1, type=int, help='File index start.')
args = parser.parse_args()

eosDir='/eos/uscms/store/user/lpcml/mandrews/IMG'
outDir='~lpcml/nobackup/mandrews' # NOTE: Space here is limited, transfer files to EOS after processing
xrootd='root://cmsxrootd.fnal.gov' # FNAL
#xrootd='root://eoscms.cern.ch' # CERN
#decays = ['QCDToGG_Pt_80_120_13TeV_TuneCUETP8M1_noPU', 'QCDToQQ_Pt_80_120_13TeV_TuneCUETP8M1_noPU']
decays = ['','QCDToQQ_Pt_80_120_13TeV_TuneCUETP8M1_noPU']


scale = [1., 1.]
chunk_size = 200
jet_shape = 125
njets = 2

# Loop over decays
for d, decay in enumerate(decays):

    if d != args.label:
        pass
        continue

    print '>> Doing decay[%d]: %s'%(d, decay)

    tfile_idxs = glob.glob('%s/%s*_IMG/*/*/output_*.root'%(eosDir,decay))
    #tfile_idxs = glob.glob('%s/%s*_IMG/*/*/output_1.root'%(eosDir,decay))
    tfile_idxs = [s.replace('.root','').split('_')[-1] for s in tfile_idxs]
    tfile_idxs = [int(i) for i in tfile_idxs]
    tfile_idxs.sort()
    #tfile_idxs = [1] # DEBUG mode: for single, local file
    print '>> File idxs:', tfile_idxs

    # Loop over root ntuples
    for n in tfile_idxs:

        if n < args.file_idx_start:
            continue

        tfile_str = glob.glob('%s/%s*_IMG/*/*/output_%d.root'%(eosDir,decay,n))
        assert len(tfile_str) == 1, "More than 1 file of same name found in different dirs: %s"%tfile_str
        tfile_str = '%s/%s'%(xrootd,tfile_str[0])
        #tfile_str = 'output_dijet.root' # DEBUG mode: for single, local file
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

        # eventId
        branches = ["eventId"]
        eventId = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.int32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],eventId.shape)

        ## lumiId
        #branches = ["lumiId"]
        #lumiId = da.concatenate([\
        #            da.from_delayed(\
        #                load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
        #                shape=(get_chunk_size(i,neff,chunk_size),),\
        #                dtype=np.int32)\
        #            for i in proc_range])
        #print " >> %s: %s"%(branches[0],lumiId.shape)

        # runId
        branches = ["runId"]
        runId = da.concatenate([\
                    da.from_delayed(\
                        load_single(tree,i,i+get_chunk_size(i,neff,chunk_size), branches),\
                        shape=(get_chunk_size(i,neff,chunk_size),),\
                        dtype=np.int32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],runId.shape)

        # ECAL
        readouts = [280,360]
        branches = ["ECAL_energy"]
        X_ECAL = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],X_ECAL.shape)

        # ECAL with resampled EE
        X_ECAL_EEup = X_ECAL.map_blocks(lambda x: block_resample_EE(x), dtype=np.float32)
        print " >> %s: %s"%('ECAL_EEup_energy',X_ECAL_EEup.shape)

        # Tracks at ECAL
        readouts = [280,360]
        branches = ["ECAL_tracksPt"]
        #branches = ["ECAL_tracksQPt"] # for Qxpt weighted 
        #branches = ["ECAL_EndtracksPt"] # for pt weighted at ECAL face
        X_TracksAtECAL = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],X_ECAL.shape)

        # HBHE upsample
        readouts = [56,72]
        branches = ["HBHE_energy"]
        upscale = 5
        X_HBHE_up = da.concatenate([\
                    da.from_delayed(\
                        load_X_upsampled(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1], upscale),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s(upsampled): %s"%(branches[0],X_HBHE_up.shape)
        
        #X_MuonsAtECAL, 
        X_ECAL_stacked = da.concatenate([X_TracksAtECAL, X_ECAL_EEup, X_HBHE_up], axis=-1)
        print " >> %s: %s"%('X_ECAL_stacked', X_ECAL_stacked.shape)

        # EB
        readouts = [170,360]
        #branches = ["HBHE_energy_EB"]
        #branches = ["TracksQPt_EB","EB_energy"]
        #branches = ["TracksPt_EB","EB_energy"]
        branches = ["EndTracksPt_EB","EB_energy"]
        #branches = ["EB_energy"]
        #branches = ["EB_energy","HBHE_energy_EB","Tracks_EB"]
        X_EB = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[0]),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],X_EB.shape)

        # EE-
        readouts = [100,100]
        branches = ["TracksPt_EEm","EEm_energy","HBHE_energy_EEm"]
        #branches = ["EndTracksPt_EEm","EEm_energy","HBHE_energy_EEm"] # for pt weighted at ECAL face
        #branches = ["TracksQPt_EEm","EEm_energy","HBHE_energy_EEm"] # for Qxpt weighted 
        #branches = ["EEm_energy","HBHE_energy_EEm","Tracks_EEm"]
        X_EEm = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],X_EEm.shape)

        # EE+
        readouts = [100,100]
        branches = ["TracksPt_EEp","EEp_energy","HBHE_energy_EEp"]
        #branches = ["EndTracksPt_EEp","EEp_energy","HBHE_energy_EEp"] # for pt weighted at ECAL face
        #branches = ["TracksQPt_EEp","EEp_energy","HBHE_energy_EEp"] # for Qxpt weighted
        #branches = ["EEp_energy","HBHE_energy_EEp","Tracks_EEp"]
        X_EEp = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],X_EEp.shape)

        # HBHE
        readouts = [56,72]
        branches = ["HBHE_energy"]
        X_HBHE = da.concatenate([\
                    da.from_delayed(\
                        load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s: %s"%(branches[0],X_HBHE.shape)

        # HBHE_EM
        #readouts = [56,72]
        #branches = ["HBHE_EMenergy"]
        #X_HBHE_EM = da.concatenate([\
        #            da.from_delayed(\
        #                load_X(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1]),\
        #                shape=(get_chunk_size(i,neff,chunk_size), readouts[0], readouts[1], len(branches)),\
        #                dtype=np.float32)\
        #            for i in proc_range])
        #print " >> %s: %s"%(branches[0],X_HBHE_EM.shape)

        # HB_EB upsample
        readouts = [34,72]
        branches = ["HBHE_energy_EB"]
        upscale = 5
        X_HBHE_EB_up = da.concatenate([\
                    da.from_delayed(\
                        load_X_upsampled(tree,i,i+get_chunk_size(i,neff,chunk_size), branches, readouts, scale[1], upscale),\
                        shape=(get_chunk_size(i,neff,chunk_size), readouts[0]*upscale, readouts[1]*upscale, len(branches)),\
                        dtype=np.float32)\
                    for i in proc_range])
        print " >> %s(upsampled): %s"%(branches[0],X_HBHE_EB_up.shape)

        X_EB = da.concatenate([X_EB, X_HBHE_EB_up], axis=-1)
        print " >> %s: %s"%('X_EB', X_EB.shape)

        # Loop over jets
        for ijet in range(njets):

            print ' >> jet index:',ijet

            # jet m0
            branches = ["jetM"]
            jetM = da.concatenate([\
                        da.from_delayed(\
                            load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                            shape=(get_chunk_size(i,neff,chunk_size),),\
                            dtype=np.float32)\
                        for i in proc_range])
            print "  >> jetM:", jetM.shape

            # jet pt
            branches = ["jetPt"]
            jetPt = da.concatenate([\
                        da.from_delayed(\
                            load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                            shape=(get_chunk_size(i,neff,chunk_size),),\
                            dtype=np.float32)\
                        for i in proc_range])
            print "  >> jetPt:", jetPt.shape

            # jet seed iphi
            branches = ["jetSeed_iphi"]
            jetSeed_iphi = da.concatenate([\
                        da.from_delayed(\
                            load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                            shape=(get_chunk_size(i,neff,chunk_size),),\
                            dtype=np.float32)\
                        for i in proc_range])
            print "  >> jetSeed_iphi:", jetSeed_iphi.shape

            # jet seed ieta
            branches = ["jetSeed_ieta"]
            jetSeed_ieta = da.concatenate([\
                        da.from_delayed(\
                            load_vector(tree,i,i+get_chunk_size(i,neff,chunk_size),branches,ijet),\
                            shape=(get_chunk_size(i,neff,chunk_size),),\
                            dtype=np.float32)\
                        for i in proc_range])
            print "  >> jetSeed_ieta:", jetSeed_ieta.shape

            # jet window
            X_jets = da.concatenate([\
                        da.from_delayed(\
                            crop_jet_block(X_ECAL_stacked[i:i+get_chunk_size(i,neff,chunk_size)],\
                                             jetSeed_iphi[i:i+get_chunk_size(i,neff,chunk_size)],\
                                             jetSeed_ieta[i:i+get_chunk_size(i,neff,chunk_size)], jet_shape),\
                            shape=(get_chunk_size(i,neff,chunk_size), jet_shape, jet_shape, X_ECAL_stacked.shape[-1]),\
                            dtype=np.float32)\
                        for i in proc_range])

            # Class label
            label = d
            print "  >> Class label:",label
            y = da.from_array(\
                    np.full(len(X_jets), label, dtype=np.float32),\
                    chunks=(get_chunk_size(i,neff,chunk_size),))
            print "  >> y shape:",y.shape

            outPath = '%s/%s_IMGjet'%(outDir, decay)
            if not os.path.isdir(outPath):
                os.makedirs(outPath)
            file_out_str = "%s/%s_IMGjet_RH%d_n%d_label%d_jet%d_%d.hdf5"%(outPath,decay,int(scale[0]),neff,label,ijet,n)
            if os.path.isfile(file_out_str):
                os.remove(file_out_str)
            print "  >> Writing to:", file_out_str
            da.to_hdf5(file_out_str, {
                                      'eventId': eventId,
                                      #'lumiId': lumiId,
                                      'runId': runId,
                                      'X_ECAL': X_ECAL,
                                      #'X_ECAL_EEup': X_ECAL_EEup,
                                      'X_ECAL_stacked': X_ECAL_stacked,
                                      'X_EB': X_EB,
                                      'X_EEm': X_EEm,
                                      'X_EEp': X_EEp,
                                      'X_HBHE': X_HBHE,
                                      #'X_HBHE_EM': X_HBHE_EM,
                                      'X_HBHE_EB_up': X_HBHE_EB_up,
                                      'jetSeed_iphi': jetSeed_iphi,
                                      'jetSeed_ieta': jetSeed_ieta,
                                      'X_jets': X_jets,
                                      'jetPt': jetPt,
                                      'jetM': jetM,
                                      '/y': y
                                      }, compression='lzf')

            print "  >> Done.\n"

