import numpy as np
from numpy.polynomial import legendre
import pandas as pd
from scipy import optimize
from astropy.stats import sigma_clip

####
# constants
comb_beginning = 950 # for calculating readout noise
comb_ending = 1000 # for calculating readout noise
default_readout_noise = 12 #readout noise of IRD detectors: ~12e- (10min exposure)
wav_boundary_yj_and_h = 1420
gain_y = 2.99
gain_h = 2.78
uncertainty_gain = 0.2
interp_non_zero_beginning = 2 #cut +2 pix from the zero flux region at the edge of the order
interp_non_zero_ending = 9 #cut -10 pix from the zero flux region at the edge of the order
continuum_non_zero_beginning = 10 #cut +10 pix from the zero flux region at the edge of the order
continuum_non_zero_ending = 60 #cut -60 pix from the zero flux region at the edge of the order
#when fitting continuum
order_fit = 23 #he polynomial order to use
max_order_fit = 50
nsigma = [2,3] #the sigma clipping threshold: tuple (low, high)
maxniter = 3 #the maximum number of iterations to do
#when making sigma clipping
sigma_lower = 1
sigma_upper = 2
maxiters = 3
####

def fit_continuum(x, y, order=6, nsigma=[0.3,3.0], maxniter=50, fitfunc='legendre'):
    """Fit the continuum using sigma clipping

    Args:
        x: the wavelengths
        y: the log-fluxes
        order: the polynomial order to use
        nsigma: the sigma clipping threshold: tuple (low, high)
        maxniter: the maximum number of iterations to do
        fitfunc: fitting function (default is legendre)

    Returns:
        The value of the continuum at the wavelengths in x

    """
    if fitfunc==None:
        A = np.vander(x - np.nanmean(x), order+1)
    m = np.ones(len(x), dtype=bool)
    for i in range(maxniter):
        if fitfunc==None:
            w = np.linalg.solve(np.dot(A[m].T, A[m]), np.dot(A[m].T, y[m]))
            mu = np.dot(A, w)
        elif fitfunc=='legendre':
            f = legendre.Legendre.fit(x[m],y[m],order)
            mu = f(x)
        resid = y - mu
        sigma = np.sqrt(np.nanmedian(resid**2))
        m_new = (resid > -nsigma[0]*sigma) & (resid < nsigma[1]*sigma)
        if m.sum() == m_new.sum():
            m = m_new
            break
        m = m_new
    return mu, m

def continuum_rsd(rsd,npix=2048,ignore_orders=None):
    """fit continuum for rsd

    Args:
        rsd: raw spectrum detector matrix
        npix: number of pixels
        ignore_order: list of orders not to be evaluated the goodness of the fitting

    Returns:
        pandas DataFrame of continuum
    """
    order_fit_itr = order_fit
    pixels = np.arange(1,npix+1)
    orders = np.arange(1,len(rsd[0])+1)
    if ignore_orders is None:
        ignore_orders = [1,len(rsd[0])-1]

    df_continuum = pd.DataFrame([],columns=pixels,index=orders)
    order = orders[0]
    while order != orders[-1]+1:
        flat_ord = rsd[:,order-1]
        cutind = np.zeros(len(flat_ord),dtype=bool)
        cutind[10:-60] = True # mask pixels at both ends of the order
        cutind[np.where(flat_ord==0)[0]] = False # mask zeros
        useind = (~np.isnan(flat_ord))
        continuum_ord, mask = fit_continuum(pixels[useind & cutind],flat_ord[useind & cutind],\
                                      order=order_fit_itr,nsigma=nsigma,maxniter=maxniter)
        resid = flat_ord[useind & cutind][mask]-continuum_ord[mask]

        if (order not in ignore_orders) and order_fit_itr<max_order_fit and 100*np.std(resid/continuum_ord[mask])>1.5:
            # change fitting parameter(s) to avoid bad fitting
            order_fit_itr += 1
            order = orders[0]
        else:
            continuum_ord_interp = np.interp(pixels,pixels[useind & cutind],continuum_ord)
            df_continuum.loc[order,pixels] = continuum_ord_interp
            order += 1
    if order_fit_itr==max_order_fit:
        print(f'WARNING: order_fit reaches the maximum value of {max_order_fit}.')
    else:
        print(f'continuum is fitted with order_fit = {order_fit_itr}.')

    return df_continuum

def continuum_oneord(wdata,flat,order):
    """fit continuum for one order

    Args:
        wdata: the wavelength calibrated target spectrum
        flat: the wavelength calibrated FLAT
        order: order number to fit

    Returns:
        spectrum and continuum of the order
    """
    wdata_tmp = wdata[(wdata['order']==order)]
    wav_tmp = wdata_tmp['wav'].values
    flux_tmp = wdata_tmp['flux'].values
    flat_tmp = flat[(flat['order']==order)]
    ffl_tmp = flat_tmp['flux'].values
    cutind = np.zeros(len(flat_tmp),dtype=bool)
    cutind[continuum_non_zero_beginning:-continuum_non_zero_ending] = True # mask pixels at both ends of the order
    cutind[np.where(ffl_tmp==0)[0]] = False # mask zeros
    useind = (~np.isnan(ffl_tmp)) & (~np.isnan(flux_tmp))
    # CHECK: fit order, clipping sigma
    continuum, mask = fit_continuum(wav_tmp[useind & cutind],ffl_tmp[useind & cutind],order=order_fit,nsigma=nsigma,maxniter=maxniter)
    continuum = np.interp(wav_tmp[useind],wav_tmp[useind & cutind],continuum)
    return wav_tmp, flux_tmp, ffl_tmp, useind, continuum

def make_blaze(wdata,flat,readout_noise,std_order=None):
    """extract blaze function for target based on FLAT

    Args:
        wdata: the wavelength calibrated target spectrum
        flat: the wavelength calibrated FLAT
        readout_noise: the readout noise of target spectrum
        std_order: Once an order number is set, the blaze functions are standardized based on that order

    Return:
        the blaze function created by scaling FLAT by a constant
    """
    def f(scale):
        res = flux_tmp[useind] - continuum*scale
        res_clipped = sigma_clip(res,sigma_lower=sigma_lower,sigma_upper=sigma_upper,maxiters=maxiters)
        std = np.std(res_clipped)
        return std
    if not std_order is None:
        wav_tmp, flux_tmp, ffl_tmp, useind, continuum = continuum_oneord(wdata,flat,std_order)
        scale = optimize.brent(f)
    else:
        scale = 1

    orders = wdata['order'].unique()
    df_continuum = pd.DataFrame([])
    for order in orders:
        wav_tmp, flux_tmp, ffl_tmp, useind, continuum = continuum_oneord(wdata,flat,order)
        order_tmp = np.ones(len(continuum))
        order_tmp[:] = order
        data = np.array([wav_tmp[useind],order_tmp,flux_tmp[useind],ffl_tmp[useind],continuum*scale]).T
        df_tmp = pd.DataFrame(data,columns=['wav','order','flux','flat','continuum'])
        df_continuum = pd.concat([df_continuum,df_tmp],ignore_index=True)
    df_continuum['nflux'] = df_continuum['flux']/df_continuum['continuum']

    if (df_continuum['wav']>wav_boundary_yj_and_h).any():
        gain = gain_h #for H band
    else:
        gain = gain_y #for Y/J band
    df_continuum['sn_ratio'] = (gain*df_continuum['flux'])/np.sqrt(gain*df_continuum['flux'] + (gain*readout_noise)**2)
    df_continuum['tmp_uncertainty'] = np.sqrt(df_continuum['flux']/gain + readout_noise**2)
    df_continuum['uncertainty'] = df_continuum['tmp_uncertainty']/df_continuum['continuum']
    return df_continuum

def normalize(df_continuum,readout_noise):
    """normalize flux after combining all orders

    Args:
        df_continuum: pandas.DataFrame that contain the blaze function
        readout_noise: the readout noise of target spectrum

    Returns:
        pandas.DataFrame of 1D normalized spectrum
    """
    ########
    if (df_continuum['wav']>wav_boundary_yj_and_h).any():
        gain = gain_h #for H band
    else:
        gain = gain_y #for Y/J band


    def sn_ratio_uncertainty_interp(df_head,df_tail,readout_noise):
        """calculate the signal-to-noise ratio and temporary uncertainty of each data point.
           temporary uncertainty is a intermediate product when calculating uncertainty.

        Args:
            df_head: the spectrum of the latter order in the overlap region.
            df_tail: the spectrum of the former order in the overlap region.
            readout_noise: the readout noise of target spectrum
        Note that the wavelength range of df_tail should contain the one of df_head.
        which means max(df_tail['wav'])>max(df_head['wav']) and min(df_tail['wav'])<min(df_head['wav'])

        Returns:
            sn_ratio: signal-to-noise ratio of each data point
            tmp_uncertainty: the uncertainty of flux before normalization
        """
        wav_tail = df_tail['wav'].values
        wav_head = df_head['wav'].values
        sn_ratio = np.zeros(len(wav_head))
        tmp_uncertainty = np.zeros(len(wav_head))
        #print('former tail: ','length: ',len(wav_tail),wav_tail[:3],wav_tail[-3:])
        #print('latter head: ','length: ',len(wav_head),wav_head[:3],wav_head[-3:])
        for i in range(len(wav_head)):
            if (wav_head[i] == wav_tail).any():
                print('Found an equal value!')
                sn_ratio[i] = df_tail['sn_ratio'].values[np.where(i == wav_tail)]
            else:
                left_wav = wav_tail[np.where(wav_tail < wav_head[i])[0][-1]]
                right_wav = wav_tail[np.where(wav_tail > wav_head[i])[0][0]]
                left_coefficient = (right_wav - wav_head[i])/(right_wav - left_wav)
                right_coefficient = (wav_head[i] - left_wav)/(right_wav - left_wav)
                left_flux = df_tail[df_tail['wav']==left_wav]['flux'].values
                right_flux = df_tail[df_tail['wav']==right_wav]['flux'].values
                if left_flux/gain + readout_noise**2 < 0.0:
                    left_uncertainty = np.nan
                else:
                    left_uncertainty = np.sqrt(left_flux/gain + readout_noise**2)
                if right_flux/gain + readout_noise**2 < 0.0:
                    right_uncertainty = np.nan
                else:
                    right_uncertainty = np.sqrt(right_flux/gain + readout_noise**2)
                #this uncertainty means the uncertainty of flux before normalization
                sn_ratio[i] = np.nan
                if np.isnan(left_uncertainty) or np.isnan(right_uncertainty):
                    tmp_uncertainty[i] = np.nan
                else:
                    tmp_uncertainty[i] = np.sqrt((left_coefficient*left_uncertainty)**2 + (right_coefficient*right_uncertainty)**2)
        return sn_ratio, tmp_uncertainty
    orders=df_continuum['order'].unique()
    # interpolate flux
    df_interp = pd.DataFrame([],columns=['wav','flux','continuum','sn_ratio','tmp_uncertainty'])
    for order in orders:
        df_form = df_continuum[df_continuum['order']==order]
        df_latt = df_continuum[df_continuum['order']==(int(order)+1)]
        #cut +2/-10 pix from the zero flux region at the edge of the order
        ind_nzflux_form = np.where(df_form.flux!=0)[0]
        ind_nzflux_latt = np.where(df_latt.flux!=0)[0]
        df_form = df_form[min(ind_nzflux_form)+interp_non_zero_beginning:max(ind_nzflux_form)-interp_non_zero_ending]
        try:
            df_latt = df_latt[min(ind_nzflux_latt)+interp_non_zero_beginning:max(ind_nzflux_latt)-interp_non_zero_ending]
        except:
            pass

        if order==min(orders):
            wav_tail = df_form['wav'].values[-1]
            wav_head = df_latt['wav'].values[0]
            df_interp = pd.concat([df_interp,df_form[df_form['wav']<wav_head][df_interp.columns]])
        elif order==max(orders):
            add_ind = (df_form['wav']>df_interp['wav'].values[-1])
            df_interp = pd.concat([df_interp,df_form[add_ind][df_interp.columns]])
            continue
        else:
            wav_tail = df_form['wav'].values[-1]
            wav_head = df_latt['wav'].values[0]
            add_ind = (df_form['wav']>df_interp['wav'].values[-1]) & (df_form['wav']<wav_head)
            df_interp = pd.concat([df_interp,df_form[add_ind][df_interp.columns]])
        df_head = df_latt[df_latt['wav']<=wav_tail]
        tmp_wav = df_form[wav_head>df_form['wav']]['wav'].values[-1]
        df_tail = df_form[tmp_wav<=df_form['wav']]
        flux_interp = np.interp(df_head['wav'].values,df_tail['wav'].values,df_tail['flux'].values)
        continuum_interp = np.interp(df_head['wav'].values,df_tail['wav'].values,df_tail['continuum'].values)
        #
        sn_ratio_interp, tmp_uncertainty = sn_ratio_uncertainty_interp(df_head,df_tail,readout_noise)
        #
        flux_sum = df_head['flux'].values+flux_interp
        continuum_sum = df_head['continuum'].values+continuum_interp
        sn_ratio_sum = sn_ratio_interp
        uncertainty_ratio_square_sum = np.sqrt(df_head['tmp_uncertainty'].values**2 + tmp_uncertainty**2)
        data = np.array([df_head['wav'].values,flux_sum,continuum_sum,sn_ratio_sum,uncertainty_ratio_square_sum]).T
        df_tmp = pd.DataFrame(data,columns=df_interp.columns)
        df_interp = pd.concat([df_interp,df_tmp])

    zeroind = df_interp['continuum']==0
    df_interp = df_interp.assign(nflux=np.nan,uncertainty=np.nan)
    df_interp.loc[~zeroind,'nflux'] = df_interp['flux'][~zeroind]/df_interp['continuum'][~zeroind]
    df_interp.loc[~zeroind,'uncertainty'] = df_interp['tmp_uncertainty'][~zeroind]/df_interp['continuum'][~zeroind]
    return df_interp

def comb_norm(wfile,flatfile,combfile=None,method='mad',blaze=True):
    """read .dat file and make 1D normalized spectrum

    Args:
        wfile: path to the wavelength calibrated 1D target spectrum
        flatfile: path to the wavelength calibrated 1D FLAT or blaze file
        combfile: path to the laser frequency laser spectrum in Y/J band corresponding to wfile.
        blaze: True when self.apext_flatfield() is used. if False, blaze function will be created from 1D FLAT

    Returns:
        pandas.DataFrame of 1D normalized spectrum
    """
    wdata = pd.read_csv(wfile,header=None,delim_whitespace=True,names=['wav','order','flux'])
    flat = pd.read_csv(flatfile,header=None,delim_whitespace=True,names=['wav','order','flux'])
    if combfile != None:
        comb = pd.read_csv(combfile,header=None,delim_whitespace=True,names=['wav','order','flux'])
        comb_tmp = comb[(comb['wav'] > comb_beginning) & (comb['wav'] < comb_ending)]
        if method=='mad':
            readout_noise = 1.4826 * np.median(np.abs(comb_tmp['flux'].values - np.median(comb_tmp['flux'].values)))
        elif method=='std':
            readout_noise = np.std(comb_tmp['flux'].values)
        print('Readout Noise is :', readout_noise)
    else:
        readout_noise = default_readout_noise
        print('Using default readout Noise :', readout_noise)
        print('readout noise of IRD detectors: ~12e- (10min exposure)')
    if blaze:
        df_continuum = blaze_to_df(wdata,flat,readout_noise)
    else:
        df_continuum = make_blaze(wdata,flat,readout_noise)
    df_interp = normalize(df_continuum,readout_noise)
    return df_continuum, df_interp

def blaze_to_df(wdata,blaze,readout_noise):
    """divide target flux by the blaze function

    Args:
        wdata: the wavelength calibrated 1D target spectrum
        blaze: blaze function

    Returns:
        normalized spectrum of each order
    """
    ### readout_noise
    blaze = blaze.assign(blaze_continuum=0)
    blaze = blaze.rename(columns={'flux':'continuum'})
    df_continuum = pd.merge(wdata,blaze,on=['wav','order'])
    df_continuum['nflux'] = df_continuum['flux']/df_continuum['continuum']
    if (df_continuum['wav']>wav_boundary_yj_and_h).any():
        gain = gain_h #for H band
    else:
        gain = gain_y #for Y/J band
    df_continuum['sn_ratio'] = (gain*df_continuum['flux'])/np.sqrt(gain*df_continuum['flux'] + (gain*readout_noise)**2)
    df_continuum['tmp_uncertainty'] = np.sqrt(df_continuum['flux']/gain + readout_noise**2)
    df_continuum['uncertainty'] = df_continuum['tmp_uncertainty']/df_continuum['continuum']
    return df_continuum

if __name__ == '__main__':
    wfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/target/w41511_m2.dat'
    flatfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/flat/wflat_m2.dat'
    df_interp = comb_norm(wfile,flatfile)
    print(df_interp)
