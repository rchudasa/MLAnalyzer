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

print(torch.__version__)
use_cuda = torch.cuda.is_available()
device = torch.device("cuda:0" if use_cuda else "cpu")
torch.backends.cudnn.benchmark = True

import argparse
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-e', '--epochs',     default=10,    type=int, help='Number of training epochs.')
parser.add_argument('-l', '--lr_init',    default=5.e-4, type=float, help='Initial learning rate.')
parser.add_argument('-b', '--resblocks',  default=3,     type=int, help='Number of residual blocks.')
parser.add_argument('-a', '--load_epoch', default=0,     type=int, help='Which epoch to start training from')
parser.add_argument('-n', '--new_lr',     default=0,     type=float, help='New learning rate when loading epoch.')
parser.add_argument('-f', '--lr_factor',  default=0.2,   type=float, help='Learning rate factor')
parser.add_argument('-p', '--patience',   default=2,     type=float, help='Learning schedule patience')
parser.add_argument('-c', '--cuda',       default=0,     type=int, help='Which gpuid to use.')
args = parser.parse_args()

lr_init = args.lr_init
lr_factor = args.lr_factor
new_lr = args.new_lr
patience = args.patience
resblocks = args.resblocks
epochs = args.epochs
load_epoch = args.load_epoch
os.environ["CUDA_VISIBLE_DEVICES"]=str(args.cuda)

#run_logger = False
run_logger = True
hcal_scale  = 1
ecal_scale  = 0.02
pt_scale    = 0.01
dz_scale    = 0.05
d0_scale    = 0.1
m0_scale    = 6
mass_bins = np.arange(100,5115+295,295)/1000. # for histogram in eval()
BATCH_SIZE = 256
#n_train = ( 502342 // BATCH_SIZE ) * BATCH_SIZE
n_train = ( 2310995 // BATCH_SIZE ) * BATCH_SIZE
n_val   = (  367944 // BATCH_SIZE ) * BATCH_SIZE
#n_train = ( 1024 // BATCH_SIZE ) * BATCH_SIZE
#n_val   = ( 1024 // BATCH_SIZE ) * BATCH_SIZE

#/storage/local/data1/gpuscratch/ddicroce/IMG/train_HToEleEle_m0p1To6_pT20To150_ctau0To3_eta0To1p4/
decay = 'HToEleEle_m0p1To6_pT20To150_ctau0To3_eta0To1p4'

expt_name = 'PTscale%.2f_ECALscale%.f_HCALscale%.f_AOD_m0o%.1f_ResNet_blocks%d_seedPos_MAEloss_epochs%d_BatchSize%d_from%d_ntrain%d_nval%d_run%d'%(pt_scale, ecal_scale, hcal_scale, m0_scale, resblocks, epochs, BATCH_SIZE, load_epoch, n_train, n_val, run)
#expt_name = 'PTscale%.f_ECALscale%.f_HCALscale%.f_AOD_m0o%.1f_ResNet_blocks%d_seedPos_LOGCOSHloss_epochs%d_from%d_ntrain%d_nval%d_run%d'%(pt_scale, ecal_scale, hcal_scale, m0_scale, resblocks, epochs, load_epoch, n_train, n_val, run)
#expt_name = 'PTscale%.f_ECALscale%.f_HCALscale%.f_AOD_m0o%.1f_ResNet_blocks%d_seedPos_HUBERloss_epochs%d_from%d_ntrain%d_nval%d_run%d'%(pt_scale, ecal_scale, hcal_scale, m0_scale, resblocks, epochs, load_epoch, n_train, n_val, run)
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
    if (true-pred).abs().cuda() < delta:
        loss = 0.5*((true-pred)**2).cuda()
    else:
        loss = delta*(pred-true).abs().cuda() - 0.5*(delta**2)
    return loss.mean()

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
        data['X_jet'][0] = pt_scale   * data['X_jet'][0] #Track pT
        data['X_jet'][1] = dz_scale   * data['X_jet'][1] #Track dZ sig
        data['X_jet'][2] = d0_scale   * data['X_jet'][2] #Track d0 sig
        data['X_jet'][3] = ecal_scale * data['X_jet'][3] #ECAL
        data['X_jet'][4] = hcal_scale * data['X_jet'][4] #HCAL        
        data['am']       = transform_y(np.float32(data['am']))
        data['iphi']     = np.float32(data['iphi'])/360.
        data['ieta']     = np.float32(data['ieta'])/140.
        data['label']    = self.label
        # Preprocessing
        # High Value Suppressuib
        data['X_jet'][1][data['X_jet'][1] < -1] = 0  #(20 cm)
        data['X_jet'][1][data['X_jet'][1] >  1] = 0  #(20 cm)
        data['X_jet'][2][data['X_jet'][2] < -1] = 0  #(10 cm)
        data['X_jet'][2][data['X_jet'][2] >  1] = 0  #(10 cm)
        # Zero-Suppression
        data['X_jet'][0][data['X_jet'][0] < 1.e-2] = 0. #(1 GeV)
        data['X_jet'][3][data['X_jet'][3] < 1.e-2] = 0. #(0.1 GeV)
        data['X_jet'][4][data['X_jet'][4] < 1.e-2] = 0. #(0.01 GeV)
        return dict(data)
    def __len__(self):
        return self.parquet.num_row_groups

logger('>> Experiment: %s'%(expt_name))

train_decays = glob.glob('/storage/local/data1/gpuscratch/ddicroce/IMG/train_%s/*.parquet*'%decay)
dset_train = ConcatDataset([ParquetDataset('%s'%d, i) for i,d in enumerate(train_decays)])
val_decays = glob.glob('/storage/local/data1/gpuscratch/ddicroce/IMG/val_%s/*.parquet*'%decay)
dset_val = ConcatDataset([ParquetDataset('%s'%d, i) for i,d in enumerate(val_decays)])

idxs_train = np.random.permutation(len(dset_train))
idxs_val   = np.random.permutation(len(dset_val))
logger('>> N samples: Train: %d + Val: %d'%(len(idxs_train), len(idxs_val)))
np.savez('MODELS/%s/idxs_train+val.npz'%(expt_name), idxs_train=idxs_train, idxs_val=idxs_val)

# Train dataset
#train_sampler = SubsetRandomSampler(idxs_train)
train_sampler = RandomSampler(dset_train, replacement=True, num_samples=n_train)
train_loader  = DataLoader(dataset=dset_train, batch_size=BATCH_SIZE, num_workers=10, pin_memory=True, sampler=train_sampler)
#train_loader  = DataLoader(dataset=dset_train, batch_size=BATCH_SIZE, num_workers=10, pin_memory=True)

# Val dataset
#val_sampler   = SubsetRandomSampler(idxs_val)
val_sampler   = RandomSampler(dset_val, replacement=True, num_samples=n_val)
val_loader    = DataLoader(dataset=dset_train, batch_size=BATCH_SIZE, num_workers=10, pin_memory=True, sampler=val_sampler)
#val_loader    = DataLoader(dataset=dset_val, batch_size=BATCH_SIZE, num_workers=10, pin_memory=True)

# Test sets
#sg_files = glob.glob('IMG/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_test/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_1*.parquet')
#sg_files = glob.glob('IMG/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_merged/HToTauTau_m3p6To4_pT30To120_ctau0To3_eta0To1p4_unbiased_1*.parquet')
#dset_sg = ConcatDataset([ParquetDataset('%s'%d, i) for i,d in enumerate(sg_files)])
#dset_sg = ParquetDataset('/storage/local/data1/gpuscratch/ddicroce/IMG/HToTauTau_m3p6To17_pT20To200_ctau0To3_eta0To1p4_noPix_noHCAL_fromNeg1GeV/HToTauTau_m3p6To17_pT20To200_ctau0To3_eta0To1p4.parquet.1', 1)
#sg_loader = DataLoader(dataset=dset_sg, batch_size=BATCH_SIZE, num_workers=10)
#logger('>> N test samples: sg: %d'%(len(dset_sg)))

#bg_files = glob.glob('IMG/HToTauTau_m3p6To15_pT20To200_ctau0To3_eta0To1p4_noPix_noHCAL/*.parquet.*')
#dset_bg = ConcatDataset([ParquetDataset('%s'%d, i) for i,d in enumerate(bg_files)])
#bg_loader = DataLoader(dataset=dset_bg, batch_size=BATCH_SIZE, num_workers=10)
#logger('>> N test samples: bg: %d'%(len(dset_bg)))

import torch_resnet_concat as networks
resnet = networks.ResNet(13, resblocks, [16, 32])
resnet.cuda()
optimizer = optim.Adam(resnet.parameters(), lr=lr_init)
lr_scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=lr_factor, patience=patience)

def do_eval(resnet, val_loader, mae_best, epoch, sample, tgt_label):
    global expt_name
    loss_ = 0.
    m_pred_, m_true_, mae_, mre_ = [], [], [], []
    iphi_, ieta_ = [], []
    now = time.time()
    ma_low = transform_y(3.6) # convert from GeV to network units
    for i, data in enumerate(val_loader):
        X, m0 = data['X_jet'].cuda(), data['am'].cuda()
        iphi, ieta = data['iphi'].cuda(), data['ieta'].cuda()
        logits = resnet([X, iphi, ieta])
        loss_ += mae_loss_wgtd(logits, m0).item()
        #loss_ += logcosh(logits, m0).item()
        #criterion = torch.nn.SmoothL1Loss()
        #loss_  += criterion(logits, m0)
        #loss_ += (logits, m0, 1).item()
        # Undo preproc on mass
        logits, m0 = inv_transform(logits), inv_transform(m0)
        mae = (logits-m0).abs()
        mre = (((logits-m0).abs())/m0)
        if i % 100 == 0:
            #logger("Validation %d/%d"%(i, len(val_loader)))
            #logger('Validation (%d/%d): m_pred: %s...'%(i, len(val_loader), str(np.squeeze(logits.tolist()[:5]))))
            #logger('Validation (%d/%d): m_true: %s...'%(i, len(val_loader), str(np.squeeze(m0.tolist()[:5]))))
            logger('Validation (%d/%d): Train loss:%f, mae:%f, mre:%f'%(i, len(val_loader), loss_/(i+1), mae.mean().item(), mre.mean().item() ))
        # Store batch metrics:
        #m_pred_.append(logits.tolist())
        #m_true_.append(m0.tolist())
        mae_.append(mae.tolist())
        mre_.append(mre.tolist())
        #iphi_.append(iphi.tolist())
        #ieta_.append(ieta.tolist())

    now = time.time() - now
    #m_true_ = np.concatenate(m_true_)
    #m_pred_ = np.concatenate(m_pred_)
    mae_    = np.concatenate(mae_)
    mre_    = np.concatenate(mre_)
    #iphi_   = np.concatenate(iphi_)
    #ieta_   = np.concatenate(ieta_)

    #logger('%d: Val m_pred: %s...'%(epoch, str(np.squeeze(m_pred_[:5]))))
    #logger('%d: Val m_true: %s...'%(epoch, str(np.squeeze(m_true_[:5]))))
    #logger('%d: Val time:%.2fs in %d steps for N=%d'%(epoch, now, len(val_loader), len(m_true_)))
    logger('%d: Val loss:%f, mae:%f, mre:%f'%(epoch, loss_/len(val_loader), np.mean(mae_), np.mean(mre_)))

    score_str = 'epoch%d_%s_mae%.4f'%(epoch, sample, np.mean(mae_))

    lr_scheduler.step(loss_/len(val_loader))
    print(optimizer.param_groups[0]['lr'])

    #if epoch == 1 or epoch == 5 or epoch == 10 or epoch == 20:
    if epoch == 100:

        if 'pseudoscalar' in sample:

            # Check 2D m_true v m_pred
            logger('%d: Val m_true vs. m_pred, [3600,12000,400] MeV:'%(epoch))
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
        logger('%d: Val m_pred, [3600,12000,400] MeV: %s'%(epoch, str(np.uint(hst))))
        mlow = hst[0]
        mrms = np.std(hst)
        logger('%d: Val m_pred, [3600,12000,400] MeV: low:%d, rms: %f'%(epoch, mlow, mrms))
        norm = 1.*len(m_pred_)/len(m0)
        plt.hist(m_true_, range=(2.4,14), bins=29, histtype='step', label=r'$\mathrm{m_{true}}$', linestyle='--', color='grey', alpha=0.6)
        plt.hist(m_pred_, range=(2.4,14), bins=29, histtype='step', label=r'$\mathrm{m_{pred}}$', linestyle='--', color='C0', alpha=0.6)
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
    epoch_string = 'MODELS/HToEleEle_m0p1To6_pT20To150_ctau0To3_eta0To1p4_PTscale0.01_ECALscale0_HCALscale1_AOD_m0o6.0_ResNet_blocks3_seedPos_MAEloss_epochs10_BatchSize256_from0_ntrain2310912_nval367872_run0/model_epoch%d'%(load_epoch) #loading validation model
    for model_name in glob.glob('%s*pkl'%(epoch_string)):
        print(model_name)
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

        if new_lr != 0:
            logger(' - OLD LR = %f'%optimizer.param_groups[0]['lr'])
            optimizer.param_groups[0]['lr'] = new_lr
            logger(' - NEW LR = %f'%optimizer.param_groups[0]['lr'])


#print_step = 10
print_step = 100
#print_step = 1000
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
        X, m0 = data['X_jet'].cuda(), data['am'].cuda()
        iphi, ieta = data['iphi'].cuda(), data['ieta'].cuda()
        optimizer.zero_grad()
        logits = resnet([X, iphi, ieta])
        loss = mae_loss_wgtd(logits, m0)
        #loss = logcosh(logits, m0)
        #loss = huber(logits, m0, 1)
        #criterion = torch.nn.SmoothL1Loss()
        #loss      = criterion(logits, m0)
        #break
        loss.backward()
        optimizer.step()
        epoch_wgt += len(m0) 
        #epoch_wgt += wgts.sum()
        n_trained += 1
        if i % print_step == 0:
            logits, m0 = inv_transform(logits), inv_transform(m0)
            mae =  (logits-m0).abs().mean()
            mre = (((logits-m0).abs())/m0).mean()
            logger('%d: (%d/%d) m_pred: %s...'%(epoch, i, len(train_loader), str(np.squeeze(logits.tolist()[:5]))))
            logger('%d: (%d/%d) m_true: %s...'%(epoch, i, len(train_loader), str(np.squeeze(m0.tolist()[:5]))))
            logger('%d: (%d/%d) Train loss:%f, mae:%f, mre:%f'%(epoch, i, len(train_loader), loss.item(), mae.item(), mre.item() ))

    now = time.time() - now
    logits, m0 = inv_transform(logits), inv_transform(m0)
    mae = (logits-m0).abs().mean()
    mre = ((logits-m0).abs()/m0).mean()
    logger('%d: Train time:%.2fs in %d steps for N:%d, wgt: %.f'%(epoch, now, len(train_loader), n_trained, epoch_wgt))
    logger('%d: Train loss:%f, mae:%f, mre:%f'%(epoch, loss.item(), mae.item(), mre.item() ))
  
    #lr_scheduler.step(loss.item())
    #logger('Adaptive LR %f'%optimizer.param_groups[0]['lr'])

    #if run_logger:
    #    loss_value = loss.item()
    #    epoch_str = 'epoch%d_mae%f'%(epoch, mae.item())
    #    filename  = 'MODELS/%s/train_model_%s.pkl'%(expt_name, epoch_str)
    #    model_dict = {'model_state_dict': resnet.state_dict(), 'optimizer_state_dict': optimizer.state_dict(), 'epoch' : epoch, 'loss': loss_value}
    #    torch.save(model_dict, filename)

    # Run Validation
    resnet.eval()
    _ = do_eval(resnet, val_loader, mae_best, epoch, 'val_pseudoscalar', 1)
    #_ = do_eval(resnet, val_loader, mae_best, epoch, 'val_tau', 0)

    #_ = do_eval(resnet, sg_loader, mae_best, epoch, 'test_pseudoscalar_M3.6To4', 1)
    #_ = do_eval(resnet, bg_loader, mae_best, epoch, 'test_pseudoscalar_M3.6To15', 1)

if run_logger:
    f.close()

