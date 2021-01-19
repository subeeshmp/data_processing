#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 10:59:54 2021

@author: subeesh
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter,filtfilt

class data_processing:

    def __init__(self, tseries, delta_t):
        self.tseries = tseries
        self.delta_t = delta_t

    # Sample Method
    def fft(self):
        pass

# butterworth filter (Added by: Vineet)
    def filter(self,fc,order=4,plot=1):
        '''
        Function to filter 1D signal
        Function picked from scipy.signal package
        
        Parameters
        ----------
        fc - cutoff frequency for the filter
        order - order of the filter, default value is 4
        plot - if 1 make plots for the original time series and filtered signal

        Returns
        -------
        fil_sig - filtered signal

        '''
        data = self.tseries
        dt   = self.delta_t
        n    = len(data)
        time = np.arange(n)*dt     # Create time axis
        n_fc = fc / (1/(2*dt))     #Normalize the frequency
        # Get the filter coefficients 
        b, a    = butter(order,n_fc,btype='low')
        fil_sig = filtfilt(b,a,data)
        # Plotting the wavelet and time series
        if plot == 1:
            fig, axs = plt.subplots(2, 1, constrained_layout=True,
                                    figsize=(8, 11))
            axs[0].plot(time, self.tseries)
            axs[0].set(title='Time series')
            axs[1].plot(time,fil_sig)
            axs[1].set(title='filtered')
            plt.show()
        return fil_sig 
