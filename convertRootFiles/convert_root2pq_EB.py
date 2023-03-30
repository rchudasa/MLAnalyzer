import ROOT
import numpy as np
import glob, os
import pyarrow as pa
import pyarrow.parquet as pq

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-i', '--infile', default=['output.root'], nargs='+', type=str, help='Input root file.')
parser.add_argument('-o', '--outdir', default='.', type=str, help='Output pq file dir.')
parser.add_argument('-d', '--decay', default='test', type=str, help='Decay name.')
parser.add_argument('-n', '--idx', default=0, type=int, help='Input root file index.')
args = parser.parse_args()

#infile_str = args.infile[0].split('/')[-4].replace('_reg','')
infile_str = args.infile[0].split('/')[-1]
print(infile_str)

rhTreeStr = args.infile
print " >> Input file[0]:",rhTreeStr[0]
rhTree = ROOT.TChain("fevt/RHTree")
for f in rhTreeStr:
  rhTree.Add(f)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> nEvts:",nEvts
outStr = '%s/%s.2016.parquet.%d'%(args.outdir, args.decay, args.idx)
print " >> Output file:",outStr

##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
#iEvtEnd   = 1000
iEvtEnd   = nEvts
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nEvts = 0
d = {} # Arrays to be written to parquet should be saved to data dict
sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    # Initialize event
    rhTree.GetEntry(iEvt)
    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    d['idx'] = [rhTree.runId, rhTree.lumiId, rhTree.eventId]
    #d['evt_wgt'] = rhTree.evt_weight

    d['hltA'] = rhTree.hltAccept

    d['mGG']  = rhTree.m0
    #if d['mGG'] < 100. or d['mGG'] > 180.: continue

    d['nRecoPho'] = rhTree.nRecoPho

    d['iphi'] = list(rhTree.SC_iphi)
    d['ieta'] = list(rhTree.SC_ieta)
    d['pho_Ecorr'] = list(rhTree.pho_ecalEPostCorr)
    d['pho_p4'] = np.transpose([rhTree.pho_E, rhTree.pho_pT, rhTree.pho_eta, rhTree.pho_phi], [1,0])
    d['pho_id'] = np.transpose([
            rhTree.pho_r9
            ,rhTree.pho_sieie
            ,rhTree.pho_phoIso
            ,rhTree.pho_chgIso
            ,rhTree.pho_chgIsoWrongVtx
            ,rhTree.pho_Eraw
            ,rhTree.pho_phiWidth
            ,rhTree.pho_etaWidth
            ,rhTree.pho_scEta
            ,rhTree.pho_sieip
            ,rhTree.pho_s4
        ], [1,0])
    d['pho_vars'] = np.transpose([
            rhTree.pho_r9
            ,rhTree.pho_HoE
            ,rhTree.pho_hasPxlSeed
            ,rhTree.pho_sieie
            ,rhTree.pho_phoIso
            ,rhTree.pho_trkIso
            ,rhTree.pho_chgIsoCorr
            ,rhTree.pho_neuIsoCorr
            ,rhTree.pho_phoIsoCorr
            ,rhTree.pho_bdt
        ], [1,0])

    #d['X'] = np.array(rhTree.SC_energy).reshape(-1,1,32,32) # nPho, nChannel, nRows, nCols
    #d['Xtz'] = np.transpose([rhTree.SC_energyT, rhTree.SC_energyZ], [1,0,2]).reshape(-1,2,32,32) # nPho, nChannel, nRows, nCols
    #d['X_aod'] = np.array(rhTree.SCaod_energy).reshape(-1,1,32,32) # nPho, nChannel, nRows, nCols
    #d['Xtz_aod'] = np.transpose([rhTree.SCaod_energyT, rhTree.SCaod_energyZ], [1,0,2]).reshape(-1,2,32,32) # nPho, nChannel, nRows, nCols
    d['X_EB'] = np.array(rhTree.EB_energy).reshape(1,170,360)

    pqdata = [pa.array([d_]) if np.isscalar(d_) or type(d_) == list else pa.array([d_.tolist()]) for d_ in d.values()]
    table = pa.Table.from_arrays(pqdata, d.keys())

    if nEvts == 0:
        writer = pq.ParquetWriter(outStr, table.schema, compression='snappy')
    writer.write_table(table)

    nEvts += 1

sw.Stop()
print " >> nDiPhoEvents:",nEvts
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
print " >> ======================================"
