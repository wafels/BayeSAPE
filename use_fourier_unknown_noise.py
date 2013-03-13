import numpy as np
import datetime
from sunpy.lightcurve import LightCurve


class BTS(input):
    """Bayesian Spectral and Parameter Estimation time-series object.
    Accepts SunPy lightcurve objects only."""
    def __init__(self, input):
        # store the light curve
        self.lightcurve = input
        
        
        # get the times from the start times
        if isinstance(self.lightcurve.data.index[0], datetime.datetime):
            self.time = (self.lightcurve.data.index[:]-self.lightcurve.data.index[0]).total_seconds()
        else:
            self.time = self.lightcurve.data.index[:]-self.lightcurve.data.index[0]

        # properties of the time-series
        self.average_cadence = self.time[:]/(1.0*(len(self.seconds)-1))
        _fftfreq=np.fft.fftfreq(len(self.seconds), self.average_cadence)
        self.fft_frequencies = _fftfreq[_fftfreq>=0.0]
        self.angular_frequency = 2*np.pi*self.frequencies

    def where_low_frequency(self, limit=0.01):
        """ Returns True for frequencies that exceed the low frequency limit"""
        return np.array(self.angular_frequency*self.time[-1]) < limit

    def get_pdf(self, use_schuster=False, frequencies=False):
        pass
        return

    def known_white_noise(self, white_noise):
        pass
        return
    def unknown_noise(d, frequency_index):
        pass