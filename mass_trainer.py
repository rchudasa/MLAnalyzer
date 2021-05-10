import numpy as np
run = 0
np.random.seed(run)
import os, glob
import time
import pyarrow as pa
import pyarrow.parquet as pq
import torch
import torch.nn.functional as F
import torch.optim as optim
import matplotlib.pyplot as plt
plt.switch_backend('agg')
#from skimage.transform import rescale
plt.rcParams["figure.figsize"] = (5,5)
from torch.utils.data import *
import pandas as pd

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-e', '--epochs', default=10, type=int, help='Number of training epochs.')
parser.add_argument('-l', '--lr_init', default=5.e-4, type=float, help='Initial learning rate.')
parser.add_argument('-b', '--resblocks', default=3, type=int, help='Number of residual blocks.')
parser.add_argument('-a', '--load_epoch', default=0, type=int, help='Which epoch to start training from')
parser.add_argument('-c', '--cuda', default=0, type=int, help='Which gpuid to use.')
args = parser.parse_args()

lr_init = args.lr_init
resblocks = args.resblocks
epochs = args.epochs
load_epoch = args.load_epoch
os.environ["CUDA_VISIBLE_DEVICES"]=str(args.cuda)

#run_logger = False
run_logger = True 
#pixel_scale = 0.2
hcal_scale  = 0.2
ecal_scale  = 0.01
pt_scale    = 0.005
m0_scale    = 14
mass_bins = np.arange(3600,12000+700,700)/1000. # for histogram in eval()
BATCH_SIZE = 64*4
#n_all = BATCH_SIZE*1594 #NO unphisical masses 
#n_val = BATCH_SIZE*500  #NO unphisical masses
n_all = BATCH_SIZE*1837 #ma>2.5 (factor 0.1)
n_val = BATCH_SIZE*500  #ma>2.5 (factor 0.1)
#n_all = BATCH_SIZE*1939 #ma>1.5 (factor 0.2)
#n_val = BATCH_SIZE*550  #ma>1.5 (factor 0.2)
#n_all = BATCH_SIZE*2074 #ma>0.5 (factor 0.3)
#n_val = BATCH_SIZE*500  #ma>0.5 (factor 0.3)
n_train = n_all - n_val

#decay = 'HToTauTau_m0p5To14_pT30To120_ctau0To3_eta0To1p4_noPix'
decay = 'HToTauTau_m2p0To14_pT30To120_ctau0To3_eta0To1p4_noPix'
#decay = 'HToTauTau_m2p5To14_pT30To120_ctau0To3_eta0To1p4_noPix'
#decay = 'HToTauTau_m3p6To14_pT30To120_ctau0To3_eta0To1p4_noPix'

#expt_name = 'PTscale%.f_ECALscale%.f_HCALscale%.f_AOD_m0o%.1f_ResNet_blocks%d_seedPos_MAEloss_epochs%d_from%d_ntrain%d_nval%d_run%d'%(pt_scale, ecal_scale, hcal_scale, m0_scale, resblocks, epochs, load_epoch, n_train, n_val, run)
expt_name = 'PTscale%.f_ECALscale%.f_HCALscale%.f_AOD_m0o%.1f_ResNet_blocks%d_seedPos_LOGCOSHloss_epochs%d_from%d_ntrain%d_nval%d_run%d'%(pt_scale, ecal_scale, hcal_scale, m0_scale, resblocks, epochs, load_epoch, n_train, n_val, run)
#expt_name = 'Pixscale%.f_PTscale%.f_ECALscale%.f_HCALscale%.f_AOD_m0o%.1f_ResNet_blocks%d_seedPos_MAEloss_epochs%d_from%d_ntrain%d_nval%d_run%d'%(pt_scale, ecal_scale, hcal_scale, m0_scale, resblocks, epochs, load_epoch, n_train, n_val, run)
#expt_name = 'EBtzo%.f_AOD_m0o%.1f_ResNet_blocks%d_seedPos_MAEloss_lr%s_epochs%d_from%d_ntrain%d_nval%d_run%d'\
#            %(eb_scale, m0_scale, resblocks, str(lr_init), epochs, load_epoch, n_train, n_val, run)

expt_name = '%s_%s'%(decay, expt_name)
if run_logger:
    if not os.path.isdir('LOGS'):
        os.makedirs('LOGS')
    f = open('LOGS/%s.log'%(expt_name), 'w')
    #for d in ['MODELS', 'METRICS','PLOTS']:
    for d in ['MODELS', 'PLOTS']:
        if not os.path.isdir('%s/%s'%(d, expt_name)):
            os.makedirs('%s/%s'%(d, expt_name))

def logger(s):
    global f, run_logger
    print(s)
    if run_logger:
        f.write('%s\n'%str(s))

def mae_loss_wgtd(pred, true, wgt=1.):
    loss = wgt*(pred-true).abs().cuda()
    return loss.mean()

# huber loss
def huber(pred, true, delta):
    loss = np.where(np.abs(true-pred) < delta , 0.5*((true-pred)**2), delta*np.abs(true - pred) - 0.5*(delta**2))
    return np.sum(loss)

# log cosh loss
def logcosh(pred, true):
    #loss = torch.mean( torch.log( torch.cosh(y - y_hat) )) 
    loss = torch.log(torch.cosh(pred - true)).cuda()
    return loss.mean()

def transform_y(y):
    return y/m0_scale

def inv_transform(y):
    return y*m0_scale

class ParquetDataset(Dataset):
    def __init__(self, filename, label):
        self.parquet = pq.ParquetFile(filename)
        #self.cols = None # read all columns
        #self.cols = ['X_jet.list.item.list.item.list.item','am','apt','iphi','ieta'] 
        self.cols = ['X_jet.list.item.list.item.list.item','am','iphi','ieta'] 
        self.label = label
    def __getitem__(self, index):
        data = self.parquet.read_row_group(index, columns=self.cols).to_pydict()
        data['X_jet'] = np.float32(data['X_jet'][0])
        data['X_jet'][-1,...] = hcal_scale * data['X_jet'][-1,...] #HCAL
        data['X_jet'][-2,...] = ecal_scale * data['X_jet'][-2,...] #ECAL
        data['X_jet'][-5,...] = pt_scale   * data['X_jet'][-5,...] #Track pT
        #data['X_jet'] = np.float32(data['X_jet'][0])/ecal_scale
        data['am'] = transform_y(np.float32(data['am']))
        #data['apt'] = np.float32(data['apt'])
        data['iphi'] = np.float32(data['iphi'])/360.
        data['ieta'] = np.float32(data['ieta'])/140.
        data['label'] = self.label
        return dict(data)
    def __len__(self):
        return self.parquet.num_row_groups

logger('>> Experiment: %s'%(expt_name))

#directory='/storage/local/data1/gpuscratch/ddicroce/IMG/%s'%decay
#files=os.listdir(directory)
#df=[pd.read_parquet(os.path.join(directory,f)) for f in files]
#df=[pd.read_parquet(os.path.join(directory,f),engine='fastparquet') for f in files]
#dset_train=pd.concat(df,ignore_index=True)
#dset_train.describe()
decays = glob.glob('/storage/local/data1/gpuscratch/ddicroce/IMG/%s/*.parquet*'%decay)
dset_train = ConcatDataset([ParquetDataset('%s'%d, i) for i,d in enumerate(decays)])

idxs = np.random.permutation(len(dset_train))
logger('>> N samples: %d'%(len(idxs)))
idxs_train = idxs[:n_train]
idxs_val = idxs[n_train:n_train+n_val]
np.savez('MODELS/%s/idxs_train+val.npz'%(expt_name), idxs_train=idxs_train, idxs_val=idxs_val)
assert len(idxs_train)+len(idxs_val) <= len(idxs), '%d vs. %d'%(len(idxs_train)+len(idxs_val), len(idxs))
# Train dataset
train_sampler = SubsetRandomSampler(idxs_train)
train_loader  = DataLoader(dataset=dset_train, batch_size=BATCH_SIZE, num_workers=10, pin_memory=True, sampler=train_sampler)
# Val dataset
val_sampler = SubsetRandomSampler(idxs_val)
val_loader  = DataLoader(dataset=dset_train, batch_size=BATCH_SIZE, num_workers=10, pin_memory=True, sampler=val_sampler)
logger('>> N samples: Train: %d + Val: %d'%(len(idxs_train), len(idxs_val)))

# Test sets
sg_files = glob.glob('IMG/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_test/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_1*.parquet')
#sg_files = glob.glob('IMG/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_merged/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_1*.parquet')
dset_sg = ConcatDataset([ParquetDataset('%s'%d, i) for i,d in enumerate(sg_files)])
#dset_sg = ParquetDataset('/storage/local/data1/gpuscratch/ddicroce/IMG/HToTauTau_m3p6To17_pT20To200_ctau0To3_eta0To1p4_noPix_noHCAL_fromNeg1GeV/HToTauTau_m3p6To17_pT20To200_ctau0To3_eta0To1p4.parquet.1', 1)
sg_loader = DataLoader(dataset=dset_sg, batch_size=BATCH_SIZE, num_workers=10)
logger('>> N test samples: sg: %d'%(len(dset_sg)))

#bg_files = glob.glob('IMG/HToTauTau_m3p6To15_pT20To200_ctau0To3_eta0To1p4_noPix_noHCAL/*.parquet.*')
#dset_bg = ConcatDataset([ParquetDataset('%s'%d, i) for i,d in enumerate(bg_files)])
#bg_loader = DataLoader(dataset=dset_bg, batch_size=BATCH_SIZE, num_workers=10)
#logger('>> N test samples: bg: %d'%(len(dset_bg)))

import torch_resnet_concat as networks
resnet = networks.ResNet(5, resblocks, [16, 32])
resnet.cuda()
optimizer = optim.Adam(resnet.parameters(), lr=lr_init)
#lr_scheduler = optim.lr_scheduler.MultiStepLR(optimizer, milestones=[10,20], gamma=0.5)
lr_scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=2)

def do_eval(resnet, val_loader, mae_best, epoch, sample, tgt_label):
    global expt_name
    loss_ = 0.
    #m_pred_, m_true_, mae_, pt_, wgts_ = [], [], [], [], []
    #m_pred_, m_true_, mae_, pt_, = [], [], [], []
    m_pred_, m_true_, mae_ = [], [], []
    iphi_, ieta_ = [], []
    label_ = []
    now = time.time()
    ma_low = transform_y(3.6) # convert from GeV to network units
    for i, data in enumerate(val_loader):
        #X, m0, pt, wgts = data['Xtz_aod'].cuda(), data['m'].cuda(), data['pt'], data['w']
        #X, m0, pt = data['X_jet'].cuda(), data['am'].cuda(), data['apt']
        X, m0 = data['X_jet'].cuda(), data['am'].cuda()
        iphi, ieta = data['iphi'].cuda(), data['ieta'].cuda()
        label = data['label'].cuda()
        #print('label ', label)
        #Mass selection
        X     = X[m0[:,0]>ma_low]
        iphi  = iphi[m0[:,0]>ma_low]
        ieta  = ieta[m0[:,0]>ma_low]
        label = label[m0[:,0]>ma_low]
        m0    = m0[m0[:,0]>ma_low]
        #pt   = pt[m0[:,0]>ma_low]
        #logits = resnet(X)
        logits = resnet([X, iphi, ieta])
        #loss_ += mae_loss_wgtd(logits, m0).item()
        loss_ += logcosh(logits, m0).item()
        # Undo preproc on mass
        logits, m0 = inv_transform(logits), inv_transform(m0)
        mae = (logits-m0).abs()
        # Store batch metrics:
        m_pred_.append(logits.tolist())
        m_true_.append(m0.tolist())
        mae_.append(mae.tolist())
        iphi_.append(iphi.tolist())
        ieta_.append(ieta.tolist())
        label_.append(label.tolist())
        #pt_.append(pt.tolist())
        #wgts_.append(wgts.tolist())

    now = time.time() - now
    label_ = np.concatenate(label_)
    m_true_ = np.concatenate(m_true_)[label_==tgt_label]
    m_pred_ = np.concatenate(m_pred_)[label_==tgt_label]
    mae_ = np.concatenate(mae_)[label_==tgt_label]
    iphi_ = np.concatenate(iphi_)[label_==tgt_label]
    ieta_ = np.concatenate(ieta_)[label_==tgt_label]
    #pt_ = np.concatenate(pt_)[label_==tgt_label]
    #wgts_ = np.concatenate(wgts_)[label_==tgt_label]

    logger('%d: Val m_pred: %s...'%(epoch, str(np.squeeze(m_pred_[:5]))))
    logger('%d: Val m_true: %s...'%(epoch, str(np.squeeze(m_true_[:5]))))
    logger('%d: Val time:%.2fs in %d steps for N=%d'%(epoch, now, len(val_loader), len(m_true_)))
    logger('%d: Val loss:%f, mae:%f'%(epoch, loss_/len(val_loader), np.mean(mae_)))

    score_str = 'epoch%d_%s_mae%.4f'%(epoch, sample, np.mean(mae_))

    if epoch == 1 or epoch == 5 or epoch == 10 or epoch == 20:

        if 'pseudoscalar' in sample:

            # Check 2D m_true v m_pred
            logger('%d: Val m_true vs. m_pred, [3600,12000,700] MeV:'%(epoch))
            sct = np.histogram2d(np.squeeze(m_true_), np.squeeze(m_pred_), bins=mass_bins)[0]
            logger(np.uint(np.fliplr(sct).T))
            # Extended version
            plt.plot(m_true_, m_pred_, ".", color='black', alpha=0.1, label='MAE = %.3f GeV'%np.mean(mae_))
            plt.xlabel(r'$\mathrm{m_{label}}$', size=16)
            plt.ylabel(r'$\mathrm{m_{pred}}$', size=16)
            plt.plot((3.6, 14), (12, 12), color='r', linestyle='--', alpha=0.5)
            plt.plot((12, 12), (2.5, 14), color='r', linestyle='--', alpha=0.5)
            plt.plot((3.6, 14), (3.6, 3.6), color='r', linestyle='--', alpha=0.5)
            plt.plot((3.6, 14), (3.6, 14), color='r', linestyle='--', alpha=0.5)
            plt.xlim(3.6, 14)
            plt.ylim(2.5, 12)
            plt.legend(loc='upper left')
            plt.savefig('PLOTS/%s/mtruevpred_%s.png'%(expt_name, score_str), bbox_inches='tight')
            plt.close()
            # Truncated version
            plt.plot(m_true_, m_pred_, ".", color='black', alpha=0.125, label='MAE = %.3f GeV'%np.mean(mae_))
            plt.xlabel(r'$\mathrm{m_{label}}$', size=16)
            plt.ylabel(r'$\mathrm{m_{pred}}$', size=16)
            plt.plot((3.6, 12), (3.6, 12), color='r', linestyle='--', alpha=0.5)
            plt.xlim(3.6, 12)
            plt.ylim(3.6, 12)
            plt.legend(loc='upper left')
            plt.savefig('PLOTS/%s/mtruevpred_%s_trunc.png'%(expt_name, score_str), bbox_inches='tight')
            plt.close()

        # Check 1D m_pred
        hst = np.histogram(np.squeeze(m_pred_), bins=mass_bins)[0]
        logger('%d: Val m_pred, [3600,12000,700] MeV: %s'%(epoch, str(np.uint(hst))))
        mlow = hst[0]
        mrms = np.std(hst)
        logger('%d: Val m_pred, [3600,12000,700] MeV: low:%d, rms: %f'%(epoch, mlow, mrms))
        norm = 1.*len(m_pred_)/len(m0)
        plt.hist(m_true_, range=(2.5,14), bins=20, histtype='step', label=r'$\mathrm{m_{true}}$', linestyle='--', color='grey', alpha=0.6)
        plt.hist(m_pred_, range=(2.5,14), bins=20, histtype='step', label=r'$\mathrm{m_{pred}}$', linestyle='--', color='C0', alpha=0.6)
        #plt.hist((m_true_ if 'pseudoscalar' in sample else np.zeros_like(m_true_)),\
        #        #range=(3.6,15), bins=20, histtype='step', label=r'$\mathrm{m_{true,w}}$', color='grey', weights=wgts_*norm)
        #        range=(3.6,17), bins=20, histtype='step', label=r'$\mathrm{m_{true,w}}$', color='grey', weights=norm)
        #plt.hist(m_pred_, range=(-1,17), bins=20, histtype='step', label=r'$\mathrm{m_{pred,w}}$', color='C0', weights=wgts_*norm)
        #plt.hist(m_pred_, range=(3.6,17), bins=20, histtype='step', label=r'$\mathrm{m_{pred,w}}$', color='C0', weights=norm)
        plt.xlim(2.5,14)
        plt.xlabel(r'$\mathrm{m}$', size=16)
        if 'pseudoscalar' in sample:
            plt.legend(loc='lower center')
        else:
            plt.legend(loc='upper right')
        plt.show()
        plt.savefig('PLOTS/%s/mpred_%s.png'%(expt_name, score_str), bbox_inches='tight')
        plt.close()

    if run_logger:
        if 'pseudoscalar' in sample and 'val' in sample:
            epoch_str = 'epoch%d_%s'%(epoch, sample)
            filename  = 'MODELS/%s/model_%s.pkl'%(expt_name, score_str.replace('_val_pseudoscalar',''))
            loss_value = loss_/len(val_loader)
            model_dict = {'model_state_dict': resnet.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 'epoch' : epoch, 'loss': loss_value}
            torch.save(model_dict, filename)

    return np.mean(mae_)

# MAIN #

if load_epoch != 0:
    epoch_string = 'MODELS/HToTauTau_m2p0To14_pT30To120_ctau0To3_eta0To1p4_noPix_PTscale0_ECALscale0_HCALscale0_AOD_m0o14.0_ResNet_blocks3_seedPos_LOGCOSHloss_epochs20_from10_ntrain342272_nval128000_run0/model_epoch%d'%(load_epoch)
    for model_name in glob.glob('%s*pkl'%(epoch_string)):
        print(model_name)
        #model_name = 'MODELS/%s/model_epoch%d.pkl'%(model_directory, load_epoch)
        logger('Loading weights from %s'%model_name)
        checkpoint = torch.load(model_name)
        resnet.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        epoch = load_epoch
        #resnet.load_state_dict(checkpoint['model_state_dict'])
        #optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        #epoch = checkpoint['epoch']
        #loss = checkpoint['loss']
        #resnet.load_weights(model_name)

print_step = 100
#print_step = 10000
mae_best = 1.
logger(">> Training <<<<<<<<")
for e in range(epochs):

    epoch = e+1+load_epoch
    epoch_wgt = 0.
    n_trained = 0
    logger('>> Epoch %d <<<<<<<<'%(epoch))

    # Run training
    #lr_scheduler.step()
    resnet.train()
    now = time.time()
    for i, data in enumerate(train_loader):
        #X, m0, wgts = data['Xtz_aod'].cuda(), data['m'].cuda(), data['w'].cuda()
        X, m0 = data['X_jet'].cuda(), data['am'].cuda()
        iphi, ieta = data['iphi'].cuda(), data['ieta'].cuda()
        optimizer.zero_grad()
        #logits = resnet(X)
        logits = resnet([X, iphi, ieta])
        #loss = mae_loss_wgtd(logits, m0, wgt=wgts)
        loss = logcosh(logits, m0)
        #loss = mae_loss_wgtd(logits, m0)
        #break
        loss.backward()
        optimizer.step()
        epoch_wgt += len(m0) 
        #epoch_wgt += wgts.sum()
        n_trained += 1
        if i % print_step == 0:
            logits, m0 = inv_transform(logits), inv_transform(m0)
            mae = (logits-m0).abs().mean()
            logger('%d: (%d/%d) m_pred: %s...'%(epoch, i, len(train_loader), str(np.squeeze(logits.tolist()[:5]))))
            logger('%d: (%d/%d) m_true: %s...'%(epoch, i, len(train_loader), str(np.squeeze(m0.tolist()[:5]))))
            logger('%d: (%d/%d) Train loss:%f, mae:%f'%(epoch, i, len(train_loader), loss.item(), mae.item()))

    now = time.time() - now
    logits, m0 = inv_transform(logits), inv_transform(m0)
    mae = (logits-m0).abs().mean()
    logger('%d: Train time:%.2fs in %d steps for N:%d, wgt: %.f'%(epoch, now, len(train_loader), n_trained, epoch_wgt))
    logger('%d: Train loss:%f, mae:%f'%(epoch, loss.item(), mae.item()))

    # Run Validation
    resnet.eval()
    _ = do_eval(resnet, val_loader, mae_best, epoch, 'val_pseudoscalar', 1)
    #_ = do_eval(resnet, val_loader, mae_best, epoch, 'val_tau', 0)

    #_ = do_eval(resnet, sg_loader, mae_best, epoch, 'test_pseudoscalar_M3.6To4', 1)
    #_ = do_eval(resnet, bg_loader, mae_best, epoch, 'test_pseudoscalar_M3.6To15', 1)

if run_logger:
    f.close()


#  File "mass_trainer.py", line 316, in <module>
#    _ = do_eval(resnet, val_loader, mae_best, epoch, 'val_pseudoscalar', 1)
#  File "mass_trainer.py", line 172, in do_eval
#    m_true_ = np.concatenate(m_true_)[label_==tgt_label]
