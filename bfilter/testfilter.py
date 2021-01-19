#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run this code to test the filter module in the data_processing code
The code imports a dummy nino sst file and applies butterworth lowpass filter.
@author: amol
Modified by vineet
"""

import numpy as np
import data_processing as dp

sst = np.loadtxt('sst_nino3.dat')  # input SST time series
sst = sst - np.mean(sst)
dt = 0.5
test = dp.data_processing(sst,dt)
test.filter(fc=1/30,order=2)
