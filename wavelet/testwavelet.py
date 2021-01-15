#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run this code to test the wavelet module in the data_processing code

The code imports a dummy nino sst file and plots the wavelet.

Created on Fri Jan 15 19:05:05 2021

@author: amol
"""

import numpy as np
import sys
import os
import inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import data_processing as dp



sst = np.loadtxt('sst_nino3.dat')  # input SST time series
sst = sst - np.mean(sst)


dt = 0.25
j1 = 7/0.25
test = dp.data_processing(sst, dt)
test.wavelet(lag1=0.725, j1=j1)
