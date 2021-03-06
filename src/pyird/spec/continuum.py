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
    useind = (~np.isnan(ffl_tmp)) & (~np.isnan(flux_tmp))
    continuum = fit_continuum(wav_tmp[useind],ffl_tmp[useind],order=35,nsigma=[3,3],maxniter=3)
    return wav_tmp, flux_tmp, useind, continuum

def make_blaze(wdata,flat,std_order=11):
    """extract blaze function for target based on FLAT

    Args:
        wdata: the wavelength calibrated target spectrum
        flat: the wavelength calibrated FLAT
        std_order: standard order of reference for scaling FLAT to blaze function

    Return:
        the blaze function created by scaling FLAT by a constant
    """
    wav_tmp, flux_tmp, useind, continuum = continuum_oneord(wdata,flat,std_order)
    def f(scale):
        res = flux_tmp[useind] - continuum*scale
        res_clipped = sigma_clip(res,sigma_lower=1,sigma_upper=2,maxiters=3)
        std = np.std(res_clipped)
        return std
    scale = optimize.brent(f)

    orders = wdata['order'].unique()
    df_continuum = pd.DataFrame([])
    for order in orders:
        wav_tmp, flux_tmp, useind, continuum = continuum_oneord(wdata,flat,order)
        order_tmp = np.ones(len(continuum))
        order_tmp[:] = order
        data = np.array([wav_tmp[useind],order_tmp,flux_tmp[useind],continuum*scale]).T
        df_tmp = pd.DataFrame(data,columns=['wav','order','flux','continuum'])
        df_continuum = pd.concat([df_continuum,df_tmp],ignore_index=True)

    return df_continuum

def normalize(df_continuum):
    """normalize flux after combining all orders

    Args:
        df_continuum: pandas.DataFrame that contain the blaze function

    Returns:
        pandas.DataFrame of 1D normalized spectrum
    """
    orders=df_continuum['order'].unique()
    # interpolate flux
    df_interp = pd.DataFrame([],columns=['wav','flux','continuum'])
    for order in orders:
        df_form = df_continuum[df_continuum['order']==order]
        df_latt = df_continuum[df_continuum['order']==(int(order)+1)]
        wav_tail = df_form['wav'].values[-1]
        try:
            wav_head = df_latt['wav'].values[0]
        except:
            wav_head = None
        if order==min(orders):
            df_interp = pd.concat([df_interp,df_form[df_form['wav']<wav_head][df_interp.columns]])
        elif order==max(orders):
            add_ind = (df_form['wav']>df_interp['wav'].values[-1])
            df_interp = pd.concat([df_interp,df_form[add_ind][df_interp.columns]])
            continue
        else:
            add_ind = (df_form['wav']>df_interp['wav'].values[-1]) & (df_form['wav']<wav_head)
            df_interp = pd.concat([df_interp,df_form[add_ind][df_interp.columns]])
        df_head = df_latt[df_latt['wav']<=wav_tail]
        df_tail = df_form[wav_head<=df_form['wav']]
        flux_interp = np.interp(df_head['wav'].values,df_tail['wav'].values,df_tail['flux'].values)
        continuum_interp = np.interp(df_head['wav'].values,df_tail['wav'].values,df_tail['continuum'].values)
        flux_sum = df_head['flux'].values+flux_interp
        continuum_sum = df_head['continuum'].values+continuum_interp
        data = np.array([df_head['wav'].values,flux_sum,continuum_sum]).T
        df_tmp = pd.DataFrame(data,columns=df_interp.columns)
        df_interp = pd.concat([df_interp,df_tmp])

    df_interp['nflux'] = df_interp['flux']/df_interp['continuum']
    return df_interp

def comb_norm(wfile,flatfile):
    """read .dat file and make 1D normalized spectrum

    Args:
        wfile: path to the wavelength calibrated 1D target spectrum
        flatfile: path to the wavelength calibrated 1D FLAT

    Returns:
        pandas.DataFrame of 1D normalized spectrum
    """
    wdata = pd.read_csv(wfile,header=None,delim_whitespace=True,names=['wav','order','flux'])
    flat = pd.read_csv(flatfile,header=None,delim_whitespace=True,names=['wav','order','flux'])
    df_continuum = make_blaze(wdata,flat)
    df_interp = normalize(df_continuum)
    return df_interp

if __name__ == '__main__':
    wfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/target/w41511_m2.dat'
    flatfile = '/Users/yuikasagi/IRD/PhDwork/pyird/data/20210317/flat/wflat_m2.dat'
    df_interp = comb_norm(wfile,flatfile)
    print(df_interp)
