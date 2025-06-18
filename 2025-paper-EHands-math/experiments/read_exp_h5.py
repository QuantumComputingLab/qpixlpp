#!/usr/bin/env python3
__author__ = "Jan Balewski"
__email__ = "janstar1122@gmail.com"

import h5py
import numpy as np

inpF='ideal_exp2x_p5.h5'
inpF='arctan_ibm_aachen_23enfb.h5'

h5f = h5py.File(inpF, 'r')

# Get arrays from file
inp_data = h5f['inp_udata'][:]
mshots = h5f['raw_mshot'][:]
rec_poly = h5f['rec_poly'][:]  # shape is (2,N)
true_poly = h5f['true_poly'][:]

N = len(inp_data)
print('M: found %d samples in file %s'%(N,inpF))
print('idx  input   mshot0  mshot1    reco +/- err    true poly')
print('-------------------------------------------------------')
for i in range(N):
    print('%2d  %6.1f  %6d  %6d    %.3f +/-%.3f  %6.3f'%(
        i, inp_data[i], mshots[i,0], mshots[i,1], 
        rec_poly[0,i], rec_poly[1,i], true_poly[i]))

h5f.close() 
