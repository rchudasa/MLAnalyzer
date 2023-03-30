import numpy as np
import ROOT
from root_numpy import tree2array
from dask.delayed import delayed
import dask.array as da

#eosDir='/eos/uscms/store/user/mba2012/IMGs/HighLumi_ROOTv2'
eosDir='/eos/uscms/store/user/mba2012/IMGs'
#decays = ["H125GGgluonfusion_Pt25_Eta23_13TeV_TuneCUETP8M1_HighLumiPileUpv3", "PromptDiPhoton_MGG80toInf_Pt25_Eta23_13TeV_TuneCUETP8M1_HighLumiPileUp"]
#decays = ["H125GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3", "PromptDiPhotonAll_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUp"]
#decays = ["H125GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3", "PromptDiPhotonAll_PtHat45_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUp"]
#decays = ['H125GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ['H100GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ['H112GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ['H137GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ['H150GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ['H175GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ['H162GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ['H087GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3']
#decays = ["PromptDiPhotonAll_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUp", "H125GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3"]
#decays = ["PromptDiPhotonAll_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpV2", "H125GGgluonfusion_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUpv3"]
decays = ['PromptDiPhotonAll_MGG80toInf_Pt25_Eta14_13TeV_TuneCUETP8M1_HighLumiPileUp_HighStats']

chunk_size = 250
scale = 100.

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

def get_likelihood(m, lhood, binsLow):
  l = 1.
  if m < binsLow[0]:
    l = lhood[0]
  elif m > binsLow[-1]:
    l = lhood[-1]
  else:
    l = lhood[m >= binsLow][-1]
  return np.float32(l)

for j,decay in enumerate(decays):

    if j == 0:
        pass
        #continue

    tfile_str = '%s/%s_FEVTDEBUG_IMG.root'%(eosDir,decay)
    #tfile_str = '%s/%s_FEVTDEBUG_nXXX_IMG.root'%(eosDir,decay)
    tfile = ROOT.TFile(tfile_str)
    tree = tfile.Get('fevt/RHTree')
    nevts = tree.GetEntries()
    neff = (nevts//1000)*1000
    #neff = 256000
    #neff = 250
    print " >> Doing decay:", decay
    print " >> Input file:", tfile_str
    print " >> Total events:", nevts
    print " >> Effective events:", neff

    # EB
    readouts = [170,360]
    branches = ["EB_energy"]
    X = da.concatenate([\
                da.from_delayed(\
                    load_X(tree,i,i+chunk_size, branches, readouts, scale),\
                    shape=(chunk_size, readouts[0], readouts[1], len(branches)),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", X.shape

    # eventId
    branches = ["eventId"]
    eventId = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.int32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", eventId.shape

    # m0
    branches = ["m0"]
    m0 = da.concatenate([\
                da.from_delayed(\
                    load_single(tree,i,i+chunk_size, branches),\
                    shape=(chunk_size,),\
                    dtype=np.float32)\
                for i in range(0,neff,chunk_size)])
    print " >> Expected shape:", m0.shape

    # EB rescaled by m0
    m0_reshape = da.reshape(m0, [-1,1,1,1])
    X_m0 = scale*X/m0_reshape 
    print " >> Expected shape:", X_m0.shape

    # Likelihood weights
    if j == 0:
      #h, bins = da.histogram(m0, bins=31, range=[80., 390.]) 
      h, bins = da.histogram(m0, bins=152, range=[82., 386.]) 
    else:
      #h, bins = da.histogram(m0, bins=8, range=[70., 150.]) 
      h, bins = da.histogram(m0, bins=38, range=[74., 150.]) 
    h = h.compute()
    h = h*np.float32(neff)/np.float32(h.sum())
    binsLow = bins[:-1]
    lhood = 1./h
    lhood = lhood/lhood.sum()
    wgt = da.from_array(np.array([get_likelihood(m, lhood, binsLow) for m in m0.compute()]), chunks=(chunk_size,))
    # If you don't want to change the bkg:
    if j == 0:
      pass
      #wgt = da.from_array(np.float32(np.ones(len(wgt))), chunks=(chunk_size,))

    # Class label
    label = j
    #label = 1
    print " >> Class label:",label
    y = da.from_array(\
            np.full(X.shape[0], label, dtype=np.float32),\
            chunks=(chunk_size,))

    file_out_str = "%s/%s_IMG_RHm0_n%dk_label%d.hdf5"%(eosDir,decay,neff//1000.,label)
    #file_out_str = "%s/%s_IMG_RH%d_n%dk_label%d_wgt10GeV.hdf5"%(eosDir,decay,int(scale),neff//1000.,label)
    #file_out_str = "%s/%s_IMG_RH%d-%d_n%dk.hdf5"%(eosDir,decay,int(scale[0]),int(scale[1]),neff//1000.)
    #file_out_str = "test.hdf5"
    print " >> Writing to:", file_out_str
    da.to_hdf5(file_out_str, {'/X': X, 'X_m0': X_m0, '/y': y, 'eventId': eventId, 'm0': m0, 'wgt': wgt}, compression='lzf')
    #da.to_hdf5(file_out_str, {'/X': X, 'X_m0': X_m0, '/y': y, 'eventId': eventId, 'm0': m0}, compression='lzf')
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y, 'eventId': eventId, 'm0': m0, 'wgt': wgt}, compression='lzf')
    #da.to_hdf5(file_out_str, {'/X': X, '/y': y, 'eventId': eventId, 'm0': m0}, compression='lzf')
    #da.to_hdf5(file_out_str, {'/y': y, 'eventId': eventId, 'm0': m0}, compression='lzf')

    print " >> Done.\n"
