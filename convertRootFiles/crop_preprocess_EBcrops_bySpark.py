from pyspark import SparkContext
from pyspark.sql import SQLContext

sc = SparkContext()
sqlContext = SQLContext(sc)

import numpy as np
from scipy.ndimage import maximum_position
from pyspark.sql import Row

s = 32
crop_size = int(s*s)
w = s//2
n_rows = 170 # n_eta
n_cols = 360 # n_phi

def crop_around_max(arr,r,c):
    #global n_rows, n_cols, w
    return np.array(arr, dtype=np.float32).reshape(n_rows,n_cols)[r-w:r+w,c-w:c+w].flatten()

def process_en(en):
    nonzero = (en > 0.)
    en[nonzero] = (np.log10(en[nonzero])+1.3)/4.
    return en

def process_t(b):
    return b/50.

def log_noise(lin_noise):
    nonzero = (lin_noise > 0.)
    lin_noise[nonzero] = np.log10(lin_noise[nonzero])
    return lin_noise

def process_evt(row):
    #global crop_size
    
    ### Get channel max ###
    arr_ref = np.array(row.EB_adc6, dtype=np.float32).reshape(n_rows,n_cols)
    r, c = maximum_position(arr_ref)
    
    ### Row object can be cast as python dict ###
    ### Note down out of range maxima ###
    row_dict = row.asDict()
    if c < w or c >= n_cols-w or r < w or r >= n_rows-w:
        evt_out = {k:np.full(crop_size, -999, dtype=np.float32).tolist() for k,arr in row_dict.iteritems()}
        evt_out['keep'] = False
        return Row(**evt_out)
    
    ### Initialize output dict as cropped input Row dict ###
    evt_out = {k:crop_around_max(arr,r,c) for k,arr in row_dict.iteritems()}
    #evt_out = {k:np.array(arr, dtype=np.float32).flatten() for k,arr in row_dict.iteritems()}
    '''
    ### Process Energy ###
    dict_en = ['EBenergy', 'EBenergyRed']
    for k in dict_en:
        evt_out[k] = process_en(evt_out[k])
    
    ### Process Time ###
    dict_t = ['EBtime', 'EBtimeRed']
    for k in dict_t:
        evt_out[k] = process_t(evt_out[k])
    '''
    ### Process Digis ###
    presample = np.mean([evt_out['EB_adc0'], evt_out['EB_adc1'], evt_out['EB_adc2']], axis=0)
    #presample = log_noise(presample)
    dict_adc = ['EB_adc%d'%sample for sample in range(10)]
    for k in dict_adc:
        evt_out[k] = process_digi(evt_out[k],presample)
    
    ### Keep event ###
    ### Pyspark only accepts list types ###
    evt_out = {k:arr.tolist() for k,arr in evt_out.iteritems()}
    evt_out['keep'] = True
    return Row(**evt_out)

# Case 1
def process_digi(adc,_):
    nonzero = (adc > 0.)
    adc[nonzero] = np.log10(adc[nonzero])-2.3
    return adc


def concat(row):
    row_dict = row.asDict()
    evt_out = [row_dict['EBenergy'], row_dict['EBtime'], \
               row_dict['EBenergyRed'], row_dict['EBtimeRed'], \
               row_dict['EB_adc0'], row_dict['EB_adc1'], row_dict['EB_adc2'], \
               row_dict['EB_adc3'], row_dict['EB_adc4'], row_dict['EB_adc5'], \
               row_dict['EB_adc6'], row_dict['EB_adc7'], row_dict['EB_adc8'], row_dict['EB_adc9']]
    return Row(features=evt_out, labels=1)

df = sqlContext \
    .read.format("org.dianahep.sparkroot") \
    .load("hdfs:/cms/bigdatasci/mandrews/SinglePhotonPt50_FEVTDEBUG_n250k_IMG.root")

#n_events = df.count()
#branch_list = df.columns
#print " >> N of events:", n_events
#print " >> Input branch list:",branch_list

df_out = df.rdd.map(process_evt).toDF()
df_out = df_out.filter(df_out.keep == True).drop(df_out.keep)
df_out = df.rdd.map(concat).toDF()

df_out.write.save("hdfs:/cms/bigdatasci/mandrews/SinglePhotonPt50_IMGCROP_n250k.parquet", format="parquet")
