import numpy as np
from root_numpy import root2array, tree2array
import h5py

def x_arr(filename, tree, branch, cols):
  arr = root2array(filename, tree, branches=[branch]) 
  print arr.dtype.names
  x = np.array([np.concatenate(x).astype(np.float32).reshape(-1,cols) for x in arr])
  x = np.squeeze(x)
  return x 

hcal_rows = [17, 20]
h = h5py.File('ECALvHCAL.hdf5','w')
for r in hcal_rows:
    x = x_arr('output_numEvent1_ieta%d.root'%r, 'fevt/RHTree', 'HBHEEnergy_EB', 72)
    print x.shape
    h.create_dataset('HBHE_%d'%r, data=x)

x = x_arr('output_numEvent1_ieta17.root', 'fevt/RHTree', 'EBenergyRed', 360)
print x.shape
h.create_dataset('EB', data=x)
