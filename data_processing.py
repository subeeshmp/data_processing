#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 10:59:54 2021

@author: subeesh
"""

class data_processing: 
 
    def __init__(self, tseries,delta_t): 
        self.tseries= tseries 
        self.delta_t=delta_t
   
    # Sample Method  
    def fft(self): 
