import pyarrow.parquet as pq
import pyarrow as pa # pip install pyarrow==0.7.1
import ROOT
import numpy as np
import glob, os

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
#parser.add_argument('-i', '--infile', default='output.root', type=str, help='Input root file.')
parser.add_argument('-i', '--infile', default=['output.root'], nargs='+', type=str, help='Input root file.')
#parser.add_argument('-i', '--infile', default=['output.root'], type=list, help='Input root file.')
parser.add_argument('-o', '--outdir', default='parquet', type=str, help='Output pq file dir.')
parser.add_argument('-d', '--decay', default='test', type=str, help='Decay name.')
parser.add_argument('-n', '--idx', default=0, type=int, help='Input root file index.')
#parser.add_argument('-w', '--wgt_file', default=None, type=str, help='Weight file.')
parser.add_argument('-w', '--wgt_files', default=None, nargs='+', type=str, help='Weight file.')
args = parser.parse_args()

def crop_EBshower(imgEB, iphi, ieta, window=32):

    # NOTE: image window here should correspond to the one used in RHAnalyzer
    off = window//2
    iphi = int(iphi)+1 # seed positioned at [15,15]
    ieta = int(ieta)+1 # seed positioned at [15,15]

    # Wrap-around on left side
    if iphi < off:
        diff = off-iphi
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,-diff:],
                                   imgEB[:,ieta-off:ieta+off,:iphi+off]), axis=-1)
    # Wrap-around on right side
    elif 360-iphi < off:
        diff = off - (360-iphi)
        img_crop = np.concatenate((imgEB[:,ieta-off:ieta+off,iphi-off:],
                                   imgEB[:,ieta-off:ieta+off,:diff]), axis=-1)
    # Nominal case
    else:
        img_crop = imgEB[:,ieta-off:ieta+off,iphi-off:iphi+off]

    return img_crop

def get_weight_2d(m0, pt, m0_edges, pt_edges, wgts):
    idx_m0 = np.argmax(m0 <= m0_edges)-1
    idx_pt = np.argmax(pt <= pt_edges)-1
    #print(idx_m0, idx_pt)
    return wgts[idx_m0, idx_pt]

hmvpts, m_edgess, pt_edgess = {}, {}, {}
if args.wgt_files is not None:
    nPasses = len(args.wgt_files)
    for p,wgt_file in enumerate(args.wgt_files):
        w = np.load(wgt_file)
        hmvpts[p], m_edgess[p], pt_edgess[p] = w['mvpt'], w['m_edges'], w['pt_edges']

rhTreeStr = args.infile
print " >> Input file:",rhTreeStr
rhTree = ROOT.TChain("fevt/RHTree")
for f in rhTreeStr:
  rhTree.Add(f)
nEvts = rhTree.GetEntries()
assert nEvts > 0
print " >> nEvts:",nEvts
outStr = '%s/%s.parquet.%d'%(args.outdir, args.decay, args.idx)
print " >> Output file:",outStr

##### EVENT SELECTION START #####

# Event range to process
iEvtStart = 0
#iEvtEnd   = 10
iEvtEnd   = nEvts
assert iEvtEnd <= nEvts
print " >> Processing entries: [",iEvtStart,"->",iEvtEnd,")"

nPhos = 0
data = {} # Arrays to be written to parquet should be saved to data dict
sw = ROOT.TStopwatch()
sw.Start()
for iEvt in range(iEvtStart,iEvtEnd):

    # Initialize event
    rhTree.GetEntry(iEvt)

    if iEvt % 10000 == 0:
        print " .. Processing entry",iEvt

    #if rhTree.m0 < 100. or rhTree.m0 > 110.:
    #  continue

    idx = [rhTree.runId, rhTree.lumiId, rhTree.eventId]

    #EB_time = np.array(rhTree.EB_time).reshape(170,360)
    #TracksAtEB_pt = np.array(rhTree.TracksPt_EB).reshape(170,360)
    #X_EB = np.stack([TracksAtEB_pt, EB_time], axis=0)
    #X_EB = np.array(rhTree.EB_energy).reshape(1,170,360)
    X_EB = np.array(rhTree.TracksPt_EB).reshape(1,170,360)

    nPhoEvt = len(rhTree.SC_mass)
    if args.wgt_files is not None:
        rands = np.random.random((nPhoEvt, nPasses))

    for i in range(nPhoEvt):

        data['idx'] = idx + [i]
        #data['m0'] = rhTree.m0

        data['m'] = rhTree.SC_mass[i]
        data['pt'] = rhTree.SC_pT[i]
        data['iphi'] = rhTree.SC_iphi[i]
        data['ieta'] = rhTree.SC_ieta[i]

        if data['ieta'] >= 170-16:
            continue

        if rhTree.pho_pT[i] > 100.: continue

        if args.wgt_files is not None:
            keepEG = True
            for p in range(nPasses):
                if rands[i,p] < get_weight_2d(data['m'], data['pt'], m_edgess[p], pt_edgess[p], hmvpts[p]):
                    keepEG = False
            if keepEG == False: continue
            #if rands[i] < get_weight_2d(data['m'], data['pt'], m_edges, pt_edges, hmvpt):
            #    continue

        data['A_p4'] = [
            rhTree.SC_E[i]
            ,rhTree.SC_pT[i]
            ,rhTree.SC_eta[i]
            ,rhTree.SC_phi[i]
            ]
        #data['A_ancestry'] = [A_pdgId[i], A_mothPdgId[i], OutPart_pdgId[i]]
        data['A_pteta'] = [rhTree.SC_pT[i], rhTree.SC_eta[i]]

        #j = 0 if i == 1 else 1
        #if not (OutPart_pdgId[j] == 22 or abs(OutPart_pdgId[i]) <= 5 or abs(OutPart_pdgId[i]) == 21):
        #    continue

        data['pho_pteta'] = [rhTree.pho_pT[i], rhTree.pho_eta[i]]
        data['pho_p4'] = [rhTree.pho_E[i], rhTree.pho_pT[i], rhTree.pho_eta[i], rhTree.pho_phi[i]]
        data['pho_id'] = [
                rhTree.pho_r9[i]
                ,rhTree.pho_sieie[i]
                ,rhTree.pho_phoIso[i]
                ,rhTree.pho_chgIso[i]
                ,rhTree.pho_chgIsoWrongVtx[i]
                ,rhTree.pho_Eraw[i]
                ,rhTree.pho_phiWidth[i]
                ,rhTree.pho_etaWidth[i]
                ,rhTree.pho_scEta[i]
                ,rhTree.pho_sieip[i]
                ,rhTree.pho_s4[i]
            ]
        data['pho_vars'] = [
                rhTree.pho_r9[i]
                ,rhTree.pho_HoE[i]
                ,rhTree.pho_hasPxlSeed[i]
                ,rhTree.pho_sieie[i]
                ,rhTree.pho_phoIso[i]
                ,rhTree.pho_trkIso[i]
                ,rhTree.pho_chgIsoCorr[i]
                ,rhTree.pho_neuIsoCorr[i]
                ,rhTree.pho_phoIsoCorr[i]
                ,rhTree.pho_bdt[i]
            ]

        data['X'] = np.array(rhTree.SC_energy[i]).reshape(1,32,32)
        sc_energyT = np.array(rhTree.SC_energyT[i]).reshape(1,32,32)
        sc_energyZ = np.array(rhTree.SC_energyZ[i]).reshape(1,32,32)
        data['Xtz'] = np.concatenate((sc_energyT, sc_energyZ), axis=0)

        #data['X_aod'] = np.array(rhTree.SCaod_energy[i]).reshape(1,32,32)
        #scaod_energyT = np.array(rhTree.SCaod_energyT[i]).reshape(1,32,32)
        #scaod_energyZ = np.array(rhTree.SCaod_energyZ[i]).reshape(1,32,32)
        #data['Xtz_aod'] = np.concatenate((scaod_energyT, scaod_energyZ), axis=0)

        #data['X_reco'] = np.array(rhTree.SCreco_energy[i]).reshape(1,32,32)
        #screco_energyT = np.array(rhTree.SCreco_energyT[i]).reshape(1,32,32)
        #screco_energyZ = np.array(rhTree.SCreco_energyZ[i]).reshape(1,32,32)
        #data['Xtz_reco'] = np.concatenate((screco_energyT, screco_energyZ), axis=0)

        sc_cms = crop_EBshower(X_EB, data['iphi'], data['ieta'])
        #if sc_cms.shape != data['Xtz'].shape:
        #if sc_cms.shape[1] != data['Xtz'].shape[1]:
        #    print(sc_cms.shape, data['ieta'], data['iphi'])
        #    print(data['Xtz'].shape)
        data['Xtzk'] = np.concatenate([sc_cms, data['Xtz']], axis=0)

        pqdata = [pa.array([d]) if np.isscalar(d) or type(d) == list else pa.array([d.tolist()]) for d in data.values()]
        table = pa.Table.from_arrays(pqdata, data.keys())

        if nPhos == 0:
            writer = pq.ParquetWriter(outStr, table.schema, compression='snappy')

        writer.write_table(table)

        nPhos += 1

writer.close()

sw.Stop()
print " >> nPhos:",nPhos
print " >> Real time:",sw.RealTime()/60.,"minutes"
print " >> CPU time: ",sw.CpuTime() /60.,"minutes"
print " >> ======================================"

pqIn = pq.ParquetFile(outStr)
print(pqIn.metadata)
print(pqIn.schema)
X = pqIn.read_row_group(0, columns=['idx.list.item','iphi','ieta','pho_p4.list.item']).to_pydict()
X = pqIn.read_row_group(1, columns=['idx.list.item','iphi','ieta','pho_p4.list.item']).to_pydict()
print(X)
#X = pqIn.read_row_group(0, columns=['X.list.item.list.item.list.item']).to_pydict()['X']
#X = pqIn.read(['X.list.item.list.item.list.item']).to_pydict()['X']
#X = np.float32(X)
