#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 10:59:54 2021

@author: subeesh
"""

import numpy as np
from wavelet.waveletFunctions import wavelet, wave_signif
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


class data_processing:

    def __init__(self, tseries, delta_t):
        self.tseries = tseries
        self.delta_t = delta_t

    # Sample Method
    def fft(self):
        pass

    def wavelet(self, pad=1, dj=0.25,
                s0=None, j1=None, mother='MORLET',
                lag1=0, sigtest=0, siglvl=0.95, plot=1):
        '''
        Function to calculate and plot wavelet power of a 1-d time-series

        The original Matlab code written by C. Torrence is
        modified to Python by Evgeniya Predybaylo, December 2014

        Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
        Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.


        Parameters
        ----------
        pad : TYPE, optional
              If set to 1, pad time series
              with enough zeroes to get N up to the next
              higher power of 2. This prevents wraparound
              from the end of the time series to the beginning,
              and also speeds up the FFT's used to do the wavelet
              transform. This will not eliminate all edge
              effects (see COI below). The default is 0.
        dj : TYPE, optional
             The spacing between discrete scales. Default is 0.25.
             A smaller # will give better scale resolution, but be
             slower to plot. The default is 0.25.
        s0 : TYPE, optional
             The smallest scale of the wavelet.
             The default is 2*DT.
        j1 : TYPE, optional
             the # of scales minus one. Scales range from
             S0 up to S0*2**(J1*DJ), to give a total of (J1+1)
             scales. Default is J1 = 7/DJ (or J1 = (LOG2(N DT/S0))/DJ).
        mother : TYPE, optional
             A string, equal to 'MORLET' or 'PAUL' or 'DOG'.
             The default is 'MORLET'.
        lag1 : LAG 1 Autocorrelation, used for SIGNIF levels.
             Default is 0.0
        sigtest :   0, 1, or 2.    If omitted, then assume 0.

             If 0 (the default), then just do a regular chi-square test,
              i.e. Eqn (18) from Torrence & Compo.
             If 1, then do a "time-average" test, i.e. Eqn (23).
              In this case, DOF should be set to NA, the number
              of local wavelet spectra that were averaged together.
              For the Global Wavelet Spectrum, this would be NA=N,
              where N is the number of points in your time series.
             If 2, then do a "scale-average" test, i.e. Eqns (25)-(28).
              In this case, DOF should be set to a
              two-element vector [S1,S2], which gives the scale
              range that was averaged together.
              e.g. if one scale-averaged scales between 2 and 8,
              then DOF=[2,8].
        plot : TYPE, optional
             If 1, then plot the wavelet. The default is 1.

        Returns
        -------
        wave : TYPE COMPLEX
            WAVE is the WAVELET transform of Y. This is a complex array
            of dimensions (N, J1+1). FLOAT(WAVE) gives the WAVELET amplitude,
            ATAN(IMAGINARY(WAVE), FLOAT(WAVE) gives the WAVELET phase.
            The WAVELET power spectrum is ABS(WAVE)**2.
            Its units are sigma**2 (the time series variance).
        period : TYPE FLOAT
            The vector of "Fourier" periods (in time units) that corresponds
            to the SCALEs.
        scale : TYPE FLOAT
            the vector of scale indices, given by S0*2**(j*DJ), j=0...J1
            where J1+1 is the total # of scales.
        coi : TYPE FLOAT
            Return the Cone-of-Influence, which is a vector
            of N points that contains the maximum period of useful information
            at that particular time. Periods greater than this are subject to
            edge effects.
        sigtest : TYPE FLOAT
            If 0 (the default), then just do a regular chi-square test,
            i.e. Eqn (18) from Torrence & Compo.

            If 1, then do a "time-average" test, i.e. Eqn (23).
            In this case, DOF should be set to NA, the number
            of local wavelet spectra that were averaged together.
            For the Global Wavelet Spectrum, this would be NA=N,
            where N is the number of points in your time series.

            If 2, then do a "scale-average" test, i.e. Eqns (25)-(28).
            In this case, DOF should be set to a
            two-element vector [S1,S2], which gives the scale
            range that was averaged together.
            e.g. if one scale-averaged scales between 2 and 8,
            then DOF=[2,8].

            Default is 0.
        siglvl: TYPE FLOAT
            Significance level to use. Default is 0.95
        '''

        # Inputs
        data = self.tseries
        dt = self.delta_t
        n = len(data)

        if not s0:
            s0 = 2*dt
        if not j1:
            j1 = 7/dj

        time = np.arange(n)*dt  # Create time axis
        data = data - np.mean(data)  # Remove mean
        data = data/np.std(data)  # Normalize the data

        variance = 1

        wave, period, scale, coi = wavelet(data, dt, pad, dj,
                                           s0, j1, mother)

        power = (np.abs(wave))**2  # compute wavelet power spectrum

        signif = wave_signif(([variance]), dt=dt, scale=scale,
                             sigtest=sigtest, siglvl=siglvl,
                             lag1=lag1, mother=mother)

        # expand signif --> (J+1)x(N) array
        sig95 = signif[:, np.newaxis].dot(np.ones(n)[np.newaxis, :])

        sig95 = power / sig95  # where ratio > 1, power is significant

        # Plotting the wavelet and time series

        if plot == 1:
            fig, axs = plt.subplots(2, 1, constrained_layout=True,
                                    figsize=(8, 11))
            axs[0].plot(time, self.tseries)
            axs[0].set(title='Time series')

            cs = axs[1].contourf(time, period, power)
            fig.colorbar(cs, ax=axs[1], shrink=0.9)

            axs[1].contour(time, period, sig95, [-99, 1], colors='k')
            axs[1].plot(time, coi, 'k')
            axs[1].set_yscale('log', base=2, subs=None)
            axs[1].set_ylim([np.min(period), np.max(period)])
            yax = plt.gca().yaxis
            yax.set_major_formatter(ticker.ScalarFormatter())
            axs[1].ticklabel_format(axis='y', style='plain')
            axs[1].set(ylabel='Period', title='Wavelet Power')

            plt.show()

        return wave, period, scale, coi, sig95
