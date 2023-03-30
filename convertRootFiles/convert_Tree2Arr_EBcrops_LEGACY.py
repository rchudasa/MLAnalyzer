import numpy as np
import ROOT
from root_numpy import tree2array

eosDir='/eos/cms/store/user/mandrews/ML/IMGs'
decays = {'pho': 'SinglePhotonPt50_FEVTDEBUG',
          'ele': 'SingleElectronPt50_FEVTDEBUG'}
tfiles = [ROOT.TFile('%s/%s_n250k_IMG_CROPS32_.root'%(eosDir,decay)) for decay in decays.itervalues()]
trees = {decay:tfile.Get('RHTree') for decay,tfile in zip(decays,tfiles)}

s = 32
crop_size = int(s*s)
pho_evts = trees['pho'].GetEntries()
ele_evts = trees['ele'].GetEntries()
print " >> Loaded photon tree: %d events"%pho_evts
print " >> Loaded electron tree: %d events"%ele_evts

def get_tree(tree, start, stop,isEle):
    y = []
    X = tree2array(tree, start=start, stop=stop) 
    if isEle:
        y = np.ones(len(X))
        print " >> Read-in electrons: [%d->%d)"%(start,start+len(X))
    else:
        y = np.zeros(len(X))
        print " >> Read-in photons: [%d->%d)"%(start,start+len(X))
    print " >> branches:",X.dtype.names
    X = np.array([np.concatenate(x).reshape(-1,crop_size) for x in X], dtype=np.float32)
    assert len(X) == len(y)
    return X, y

# Serve training data by chunks
# Iterates over n_chunks in data
def train_generator(n_chunks_):
    global pho_tree, ele_tree, chunk_size, s
    i = 0
    while i < n_chunks_:
        print ' >> Generating chunk: %d/%d'%(i,n_chunks_-1)
        # Get photons
        X_pho, y_pho = get_tree(trees['pho'], start=i*chunk_size, stop=(i+1)*chunk_size, isEle=False) 
        # Get electrons
        X_ele, y_ele = get_tree(trees['ele'], start=i*chunk_size, stop=(i+1)*chunk_size, isEle=True) 
        # Concatenate and reshape
        X = np.concatenate((X_pho,X_ele))
        X = X.reshape((-1,X.shape[1],s,s))
        X = np.swapaxes(X,1,2)
        X = np.swapaxes(X,2,3)
        y = np.concatenate((y_pho,y_ele))
        assert len(X) == len(y)
        yield X, y 
        i += 1

# Serve evaluation set
# Not an actual generator--serves a fixed set
def eval_generator(offset):
    global trees, chunk_size, s
    i = offset
    print ' >> Generating evalutation set:'
    # Get photons
    X_pho, y_pho = get_tree(trees['pho'], start=i*chunk_size, stop=(i+1)*chunk_size, isEle=False) 
    # Get electrons
    X_ele, y_ele = get_tree(trees['ele'], start=i*chunk_size, stop=(i+1)*chunk_size, isEle=True) 
    # Concatenate and reshape
    X = np.concatenate((X_pho,X_ele))
    X = X.reshape((-1,X.shape[1],s,s))
    X = np.swapaxes(X,1,2)
    X = np.swapaxes(X,2,3)
    y = np.concatenate((y_pho,y_ele))
    assert len(X) == len(y)
    return X, y 

#chunk_size = 25000
chunk_size = 2000
n_chunks = int(np.ceil(float(np.amax([pho_evts, ele_evts]))/float(chunk_size)))
print " >> chunk size per file:",chunk_size
print " >> total chunks:",n_chunks
epochs = 2
min_loss = 1.
for e in range(epochs):
    print " >> Epoch %d"%(e)
    for X_train, y_train in train_generator(n_chunks-2):
        print " >> X_train.shape",X_train.shape
        #model.fit(X_train, y_train, batch_size=128, nb_epoch=1, callbacks=[lr])
    # Evaluate on CV set
    X_cv, y_cv = eval_generator(n_chunks-2)
    print " >> X_cv.shape",X_cv.shape
    #cv_score = model.evaluate(X_cv, y_cv)
    #if val_loss < min_loss:
        #min_loss == val_loss
        #save model

X_test, y_cv = eval_generator(n_chunks-1)
print " >> X_test.shape",X_test.shape
