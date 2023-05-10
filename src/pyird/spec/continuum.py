import numpy as np
import pandas as pd
from scipy import optimize
from astropy.stats import sigma_clip

def fit_continuum(x, y, order=6, nsigma=[0.3,3.0], maxniter=50):
    """Fit the continuum using sigma clipping

    Args:
        x: the wavelengths
        y: the log-fluxes
        order: the polynomial order to use
        nsigma: the sigma clipping threshold: tuple (low, high)
        maxniter: the maximum number of iterations to do

    Returns:
        The value of the continuum at the wavelengths in x

    """
    A = np.vander(x - np.nanmean(x), order+1)
    m = np.ones(len(x), dtype=bool)
    for i in range(maxniter):
        w = np.linalg.solve(np.dot(A[m].T, A[m]), np.dot(A[m].T, y[m]))
        mu = np.dot(A, w)
        resid = y - mu
        sigma = np.sqrt(np.nanmedian(resid**2))
        m_new = (resid > -nsigma[0]*sigma) & (resid < nsigma[1]*sigma)
        if m.sum() == m_new.sum():
            m = m_new
            break
        m = m_new
    return mu

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
    cutind[10:-60] = True # mask pixels at both ends of the order
    useind = (~np.isnan(ffl_tmp)) & (~np.isnan(flux_tmp))
    # CHECK: fit order, clipping sigma
    continuum = fit_continuum(wav_tmp[useind & cutind],ffl_tmp[useind & cutind],order=23,nsigma=[2,3],maxniter=3)
    continuum = np.interp(wav_tmp[useind],wav_tmp[useind & cutind],continuum)
    return wav_tmp, flux_tmp, ffl_tmp, useind, continuum

def make_blaze(wdata,flat,RON,std_order=None):
    """extract blaze function for target based on FLAT

    Args:
        wdata: the wavelength calibrated target spectrum
        flat: the wavelength calibrated FLAT
        std_order: Once an order number is set, the blaze functions are standardized based on that order

    Return:
        the blaze function created by scaling FLAT by a constant
    """
    def f(scale):
        res = flux_tmp[useind] - continuum*scale
        res_clipped = sigma_clip(res,sigma_lower=1,sigma_upper=2,maxiters=3)
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
    ###############
    if (wdata['wav']>1420).any():
        gain = 2.78 #for H band
    else:
        gain = 2.99 #for Y/J band
    uncertainty_gain = 0.2
    #RON = 12 #readout noise of IRD detectors: ~12e- (10min exposure)
    df_continuum['SNratio'] = (gain*df_continuum['flux'])/np.sqrt(gain*df_continuum['flux'] + (gain*RON)**2)
    df_continuum['tmp_uncertainty'] = np.sqrt(df_continuum['flux']/gain + RON**2)
    df_continuum['uncertainty'] = df_continuum['tmp_uncertainty']/df_continuum['continuum']
    ###############
    return df_continuum

def normalize(df_continuum,RON):
    """normalize flux after combining all orders

    Args:
        df_continuum: pandas.DataFrame that contain the blaze function

    Returns:
        pandas.DataFrame of 1D normalized spectrum
    """
    ########
    if (df_continuum['wav']>1420).any():
        gain = 2.78 #for H band
    else:
        gain = 2.99 #for Y/J band
    uncertainty_gain = 0.2

    #RON = 12 #readout noise of IRD detectors: ~12e- (10min exposure)

    def SNratio_uncertainty_interp(df_head,df_tail,RON):
        wav_tail = df_tail['wav'].values
        wav_head = df_head['wav'].values
        SN = np.zeros(len(wav_head))
        tmp_uncertainty = np.zeros(len(wav_head))
        #print('former tail: ','length: ',len(wav_tail),wav_tail[:3],wav_tail[-3:])
        #print('latter head: ','length: ',len(wav_head),wav_head[:3],wav_head[-3:])
        for i in range(len(wav_head)):
            if (wav_head[i] == wav_tail).any():
                print('Found an equal value!!!')
                SN[i] = df_tail['SNratio'].values[np.where(i == wav_tail)]
            else:
                left_wav = wav_tail[np.where(wav_tail < wav_head[i])[0][-1]]
                right_wav = wav_tail[np.where(wav_tail > wav_head[i])[0][0]]
                left_coefficient = (right_wav - wav_head[i])/(right_wav - left_wav)
                right_coefficient = (wav_head[i] - left_wav)/(right_wav - left_wav)
                left_flux = df_tail[df_tail['wav']==left_wav]['flux'].values
                right_flux = df_tail[df_tail['wav']==right_wav]['flux'].values
                if left_flux/gain + RON**2 < 0.0:
                    left_uncertainty = np.nan
                else:
                    left_uncertainty = np.sqrt(left_flux/gain + RON**2)
                if right_flux/gain + RON**2 < 0.0:
                    right_uncertainty = np.nan
                else:
                    right_uncertainty = np.sqrt(right_flux/gain + RON**2)
                #this uncertainty means the uncertainty of flux before normalization
                SN[i] = np.nan
                if np.isnan(left_uncertainty) or np.isnan(right_uncertainty):
                    tmp_uncertainty[i] = np.nan
                else:
                    tmp_uncertainty[i] = np.sqrt((left_coefficient*left_uncertainty)**2 + (right_coefficient*right_uncertainty)**2)
        return SN, tmp_uncertainty
    ######
    orders=df_continuum['order'].unique()
    # interpolate flux
    df_interp = pd.DataFrame([],columns=['wav','flux','continuum','SNratio','tmp_uncertainty'])
    for order in orders:
        df_form = df_continuum[df_continuum['order']==order]
        df_latt = df_continuum[df_continuum['order']==(int(order)+1)]
        #cut +2/-10 pix from the zero flux region at the edge of the order
        ind_nzflux_form = np.where(df_form.flux!=0)[0]
        ind_nzflux_latt = np.where(df_latt.flux!=0)[0]
        df_form = df_form[min(ind_nzflux_form)+2:max(ind_nzflux_form)-9]
        try:
            df_latt = df_latt[min(ind_nzflux_latt)+2:max(ind_nzflux_latt)-9]
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
        SNratio_interp, tmp_uncertainty = SNratio_uncertainty_interp(df_head,df_tail,RON)
        #
        flux_sum = df_head['flux'].values+flux_interp
        continuum_sum = df_head['continuum'].values+continuum_interp
        SNratio_sum = SNratio_interp
        uncertainty_ratio_square_sum = np.sqrt(df_head['tmp_uncertainty'].values**2 + tmp_uncertainty**2)
        data = np.array([df_head['wav'].values,flux_sum,continuum_sum,SNratio_sum,uncertainty_ratio_square_sum]).T
        df_tmp = pd.DataFrame(data,columns=df_interp.columns)
        df_interp = pd.concat([df_interp,df_tmp])

    df_interp['nflux'] = df_interp['flux']/df_interp['continuum']
    ###
    df_interp['uncertainty'] = df_interp['tmp_uncertainty']/df_interp['continuum']
    ###
    return df_interp

def comb_norm(wfile,flatfile,LFC_path):
    """read .dat file and make 1D normalized spectrum

    Args:
        wfile: path to the wavelength calibrated 1D target spectrum
        flatfile: path to the wavelength calibrated 1D FLAT

    Returns:
        pandas.DataFrame of 1D normalized spectrum
    """
    wdata = pd.read_csv(wfile,header=None,delim_whitespace=True,names=['wav','order','flux'])
    flat = pd.read_csv(flatfile,header=None,delim_whitespace=True,names=['wav','order','flux'])
    LFC_wdat = pd.read_csv(LFC_path,header=None,delim_whitespace=True,names=['wav','order','flux'])
    df_LFC_tmp = LFC_wdat[(LFC_wdat['wav'] > 950) & (LFC_wdat['wav'] < 1000)]
    RON = 1.4826 * np.median(np.abs(df_LFC_tmp['flux'].values - np.median(df_LFC_tmp['flux'].values)))
    print('RON is :', RON)
    df_continuum = make_blaze(wdata,flat,RON)
    df_interp = normalize(df_continuum,RON)
    return df_continuum, df_interp

if __name__ == '__main__':
    wfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/target/w41511_m2.dat'
    flatfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/flat/wflat_m2.dat'
    #df_interp = comb_norm(wfile,flatfile)
    #print(df_interp)
