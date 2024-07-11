import numpy as np
import pandas as pd

__all__ = ["FluxUncertainty"]

class FluxUncertainty():
    """Calculating Uncertainty Class"""
    def __init__(self, 
                wav_boundary_yj_and_h = 1420, 
                gain_y = 2.99, 
                gain_h = 2.78, 
                combfile = None,
                comb_readout_noise_wav = [950, 1000], 
                calc_method_readout_noise = 'mad'
                ):
        """initialization of FluxUncertainty

        Args:
            wav_boundary_yj_and_h: boundary of wavelength between Y/J and H band
            gain_y: gain of IRD Y/J band detector
            gain_h: gain of IRD H band detector
            combfile: file of LFC data (mmf1)
            comb_readout_noise_wav: wavelength range [low, upp] used for calculation of readout noise
            calc_method_readout_noise: method for calculating readout noise. 'mad' is for Median Absolute Deviation, 'std' is for Standard Deviation

        """

        self.wav_boundary_yj_and_h = wav_boundary_yj_and_h 
        self.gain_y = gain_y 
        self.gain_h = gain_h 

        self.combfile = combfile
        self.comb_readout_noise_wav = comb_readout_noise_wav
        self.method = calc_method_readout_noise

        self.default_readout_noise = 12  #readout noise of IRD detectors: ~12e- (10min exposure)
        self.scale_normalized_mad = 1.4826

        self.determine_readout_noise()

    def determine_gain(self, df_continuum):
        """determine the gain to use based on wavelength

        Args:
            df_continuum: DataFrame that contains wavelength data

        Returns:
            gain of the corresponding band

        """
        if (df_continuum['wav'] > self.wav_boundary_yj_and_h).any():
            self.gain = self.gain_h #for H band
        else:
            self.gain = self.gain_y #for Y/J band

    def determine_readout_noise(self):
        """determine the readout noise 

        Returns:
            if combfile is defined, calculated real value. if not, use default value.
        """
        combfile = self.combfile
        comb_beginning, comb_ending = self.comb_readout_noise_wav

        if combfile != None:
            read_args = {'header': None, 'sep': r'\s+', 'names':['wav','order','flux']}
            comb = pd.read_csv(combfile, **read_args) 
            comb_tmp = comb[(comb['wav'] > comb_beginning) & (comb['wav'] < comb_ending)]
            if self.method=='mad':
                readout_noise = self.scale_normalized_mad * np.median(np.abs(comb_tmp['flux'].values - np.median(comb_tmp['flux'].values)))
            elif self.method=='std':
                readout_noise = np.std(comb_tmp['flux'].values)
            print('Readout Noise is :', readout_noise)
        else:
            readout_noise = self.default_readout_noise
            print('Using default readout Noise :', readout_noise)
            print('readout noise of IRD detectors: ~12e- (10min exposure)')

        self.readout_noise = readout_noise
    
    def calc_uncertainty(self, df_continuum):
        """calculation of uncertainty

        Args:
            df_continuum: DataFrame that contains flux, continuum, and tmp_uncertainty

        Returns:
            add df to columns of sn_ratio and uncertainty
        """
        self.determine_gain(df_continuum)
        
        df_continuum['sn_ratio'] = (self.gain*df_continuum['flux']) \
                                    /np.sqrt(self.gain*df_continuum['flux'] + (self.gain*self.readout_noise)**2)
        df_continuum['tmp_uncertainty'] = np.sqrt(df_continuum['flux']/self.gain + self.readout_noise**2)
        df_continuum['uncertainty'] = df_continuum['tmp_uncertainty']/df_continuum['continuum']
        return df_continuum


    def calc_uncertainty_overlap_region(self, df_head, df_tail):
        """calculate the signal-to-noise ratio and temporary uncertainty of each data point.
        
        Notes:
            the wavelength range of df_tail should contain the one of df_head.
            i.e., max(df_tail['wav'])>max(df_head['wav']) and min(df_tail['wav'])<min(df_head['wav'])
            See details in the master thesis by Ziying Gu 

        Args:
            df_head: the spectrum of the latter order in the overlap region.
            df_tail: the spectrum of the former order in the overlap region.

        Returns:
            sn_ratio: signal-to-noise ratio of each data point
            tmp_uncertainty: the uncertainty of flux before normalization 
                            (this is a intermediate product for calculating uncertainty)

        """
        if not hasattr(self, 'gain'):
            self.determine_gain(df_head)

        wav_head = df_head['wav'].values
        wav_tail = df_tail['wav'].values
        flux_tail = df_tail['flux'].values

        sn_ratio = np.zeros(len(wav_head))
        tmp_uncertainty = np.zeros(len(wav_head))

        for i in range(len(wav_head)):
            if (wav_head[i] == wav_tail).any():
                print('Found an equal value!')
                sn_ratio[i] = df_tail['sn_ratio'].values[wav_tail == wav_head[i]] 
            else:
                sn_ratio[i] = np.nan

                left_flux, right_flux, left_scale, right_scale = self.split_tail_by_wavhead(wav_tail, wav_head[i], flux_tail)
                tmp_uncertainty[i] = self.calculate_tmp_uncertainty(left_flux, right_flux, left_scale, right_scale)

        return sn_ratio, tmp_uncertainty

    def split_tail_by_wavhead(self, wav_tail, wav_head_i, flux_tail):
        """split wavelength and flux in tail by wavlength in head

        Args:
            wav_tail: wavelengths in tail (overlapping region of the former order)
            wav_head_i: wavelength in head (overlapping region of the latter order)
            flux_tail: flux in tail

        Returns:
            flux and contributions of the adjacent data points for interpolating (left and right)

        """
        left_wav = wav_tail[np.where(wav_tail < wav_head_i)[0][-1]]
        right_wav = wav_tail[np.where(wav_tail > wav_head_i)[0][0]]
        left_scale = (right_wav - wav_head_i) / (right_wav - left_wav)
        right_scale = (wav_head_i - left_wav) / (right_wav - left_wav)
        left_flux = flux_tail[wav_tail == left_wav]
        right_flux = flux_tail[wav_tail == right_wav]
        return left_flux, right_flux, left_scale, right_scale

    def calculate_tmp_uncertainty(self, left_flux, right_flux, left_scale, right_scale):
        """calculate temporal uncertainty

        Args:
            left_flux: flux at a point to the left (shorter wavelength) of the point after interpolation
            right_flux: flux at a point to the right (longer wavelength) of the point after interpolation
            left_scale: weight for left uncertainty
            right_scale: weight for right uncertainty

        Returns:
            temporary error for interpolated points

        """
        left_uncertainty = np.sqrt(left_flux / self.gain + self.readout_noise**2) \
                            if left_flux / self.gain + self.readout_noise**2 >= 0 else np.nan
        right_uncertainty = np.sqrt(right_flux / self.gain + self.readout_noise**2) \
                            if right_flux / self.gain + self.readout_noise**2 >= 0 else np.nan
  
        #this uncertainty means the uncertainty of flux before normalization
        if np.isnan(left_uncertainty) or np.isnan(right_uncertainty):
            tmp_uncertainty = np.nan
        else:
            tmp_uncertainty = np.sqrt((left_scale*left_uncertainty)**2 + (right_scale*right_uncertainty)**2)
        return tmp_uncertainty
