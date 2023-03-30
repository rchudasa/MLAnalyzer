import numpy as np
#import ROOT
#from scipy.sparse import csr_matrix
#from root_numpy import root2array, tree2array
#from root_pandas import read_root
#import tables

eosDir='/eos/cms/store/user/mandrews/ML/IMGs'

#decay = "SingleElectronPt50_FEVTDEBUG"
decay = "SinglePhotonPt50_FEVTDEBUG"
#decay = "Photon"

in_dir = 'IMGs'
s = 32
crop_size = int(s*s)

def x_generator():
    i = 0
    while i < n_chunks:
        print ' >> Chunk',i
        # Get photon
        file_in_str ='%s/%s_n250k_IMG_CROPS32_%d.npz'%(in_dir,'SinglePhotonPt50_FEVTDEBUG',i) 
        print " >> Reading",file_in_str
        with np.load(file_in_str) as file_in:
            X_pho = file_in['x']
        # Get electron
        file_in_str = '%s/%s_n250k_IMG_CROPS32_%d.npz'%(in_dir,'SingleElectronPt50_FEVTDEBUG',i) 
        print " >> Reading",file_in_str
        with np.load(file_in_str) as file_in:
            X_ele = file_in['x']
        X = np.concatenate((X_pho,X_ele))
        assert len(X) == chunk_size, "len(X) != 50k, instead %d"%len(X) 
        X = X.reshape((-1,X.shape[1],s,s))
        X = np.swapaxes(X,1,2)
        X = np.swapaxes(X,2,3)
        y = np.concatenate((np.zeros(len(X_pho),dtype=np.float32), np.ones(len(X_ele),dtype=np.float32))) 
        assert len(y) == chunk_size
        yield X, y 

def return_cv(chunk):
    print ' >> Getting CV'
    # Get photon
    file_in_str ='%s/%s_n250k_IMG_CROPS32_%d.npz'%(in_dir,'SinglePhotonPt50_FEVTDEBUG',chunk) 
    print " >> Reading",file_in_str
    with np.load(file_in_str) as file_in:
        X_pho = file_in['x']
    # Get electron
    file_in_str = '%s/%s_n250k_IMG_CROPS32_%d.npz'%(in_dir,'SingleElectronPt50_FEVTDEBUG',chunk) 
    print " >> Reading",file_in_str
    with np.load(file_in_str) as file_in:
        X_ele = file_in['x']
    X = np.concatenate((X_pho,X_ele))
    assert len(X) == chunk_size, "len(X) != 50k, instead %d"%len(X) 
    X = X.reshape((-1,X.shape[1],s,s))
    X = np.swapaxes(X,1,2)
    X = np.swapaxes(X,2,3)
    y = np.concatenate((np.zeros(len(X_pho),dtype=np.float32), np.ones(len(X_ele),dtype=np.float32))) 
    assert len(y) == chunk_size
    yield X, y 
    i += 1

#n_events_tot = 225000.
#n_events_file = 25000.
#n_chunks = int(np.ceil(n_events_tot/n_events_file))
chunk_size = 50000
n_chunks = 8
print n_chunks
epochs = 5
min_loss = 1.
for e in range(epochs):
    print " >> Epoch %d"%(e)
    for X_train, y_train in x_generator():
        pass
        print " >> X_train.shape",X_train.shape
        #model.fit(X_train, y_train, batch_size=128, nb_epoch=1, callbacks=[lr])
    # Evaluate on CV set
    #X_cv, y_cv = return_cv(8)
    #cv_score = model.evaluate(X_cv, y_cv)
    if val_loss < min_loss:
        min_loss == val_loss
        #save model
