import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

class data_processing:

    def __init__(self, tseries, delta_t):
        self.tseries = tseries
        self.delta_t = delta_t
    def npfft(self):
        '''Return the normalised amplitude of fft and frequency using fast FFT in Numpy'''
        N= self.tseries.size
        fft = np.fft.fft(self.tseries)
        ff=np.fft.fftfreq(N, d=self.delta_t)
        f=ff[:N // 2]
        abs_fft=np.abs(fft)[:N//2] * 1 / N
        return abs_fft,f
    def detrend_scipy(self):
        tseries_dtrend = signal.detrend(self.tseries)
        return tseries_dtrend
    def tapering_window(self,method="Hanning"):
        """
        Returns Tapered 1D time series.
        Input:
            1) tseries : 1D variable containing time series data.
        2) method  : Use any one of the teppering window given below.
                    "Hanning" (default)
                    "Hamming"
                    "Boxcar"
                    "Tukey"
                    "Blackman-Harris"
                    "Kaiser".
                    Output:
                        1) tseries_tapered    : Time series after applaying tapering window
    Example:
        tapering_window(time_series)
        tapering_window(time_series,"Tukey")
        """
        length=len(self.tseries)
        if method=="Hanning":
            # The Hanning window is a taper formed by using a weighted cosine.
            # w(n)=  0.5-0.5*cos(2*pi*n/(M-1)))
            #  with      0<=M<=1
            wind_hanning=np.hanning(length)
            tseries_tapered=self.tseries*wind_hanning
            print("Hanning window is applied.....")
        elif method=="Hamming":
                # The Hanning window is a taper formed by using a weighted cosine.
                # w(n)=  0.54-0.46*cos(2*pi*n/(M-1)))
                #  with      0<=M<=1
                wind_hamming=np.hamming(length)
                tseries_tapered=self.tseries*wind_hamming
                print("Hamming window is applied.....")
        elif method=="Boxcar":
            wind_boxcar=signal.boxcar(length)
            tseries_tapered=self.tseries*wind_boxcar
            print("Boxcar window is applied.....")
        elif method=="Tukey":
            # Also called "Planck-taper" window is a bump function
            wind_tukey = signal.tukey(length)
            tseries_tapered=self.tseries*wind_tukey
            print("Tukey window is applied.....")
        elif method=="Blackman-Harris":
                # A generalization of the Hamming family
                # w(n)=a0-a1cos(2*pi*n/N)+a2cos(4*pi*n/N)-a3cos(6*pi*n/N)
                wind_buckman=signal.blackmanharris(length)
                tseries_tapered=self.tseries*wind_buckman
                print("Blackman-Harris window is applied.....")
        elif method=="Kaiser":
      # The Kaiser window is a taper formed by using a Bessel function.
      # w(n)=I0(beta*sqrt(1=(4*n^2/(M-1)^2))/(I0*beta)
      # with
      #   -(M-1)/2<=n<=(M-1)/2
      # The Kaiser can approximate many other windows by varying
      # the beta parameter.
           beta=14.0
           wind_kaiser=np.kaiser(length,beta)
           tseries_tapered=self.tseries*wind_kaiser
           print("kaiser window is applied.....")
           return tseries_tapered

    
    def periodogram(self,tap_method="Hanning",deterend="cubic"):
        ''' Function for generate a periodogram 
        1) tseries : 1D variable containing time series data.
        2) Time axis
        3)tap_method  : Use any one of the teppering window given below.
                    "Hanning" (default)
                    "Hamming"
                    "Boxcar"
                    "Tukey"
                    "Blackman-Harris"
                    "Kaiser"
                    "none"
        4)deterend= Detrending method
                   "cubic"
                   "none"
        5) test= just for testing no input data required
              "yes"
               "no"
               Output:
        1) powerspectral density
        3) frequency
        2) plot
    
        Example:
        periodogram_main(tseries,t,tap_method="Hanning",deterend="cubic",test='no')   
         '''
        if deterend=="cubic":
            tseries_dtrend=data_processing.detrend_scipy(self.tseries)
        elif deterend=='none':
            tseries_dtrend=self.tseries
        if tap_method=="none":        
            tseries_tapered=tseries_dtrend
        else:
            tseries_tapered=data_processing.tapering_window(self.tseries,tap_method)
        [abs_fft,f]=data_processing.npfft(tseries_tapered,self.delta_t)
        plt.subplot(1,2,1)
        plt.plot(self.delta_t,self.tseries,'k-')
        plt.xlabel('time')
        plt.ylabel('amplitude') 
        plt.subplot(1,2,2)
        plt.plot(f,abs_fft,'k-')
        plt.xlabel('frequency')
        plt.ylabel('psd')
        plt.show() 
        return abs_fft,f       
print('check:ok')    
    
