import numpy as np
from numpy.polynomial import legendre
import pandas as pd
import warnings
    
__all__ = ["ContinuumFit"]

class ContinuumFit():
    """Continuum Fitting Class"""

    def __init__(self, 
                continuum_nonzero_ind = [10, 60], 
                base_order_fit = 23, 
                max_order_fit = 50, 
                nsigma_continuumfit = [2, 3], 
                maxiter_continuumfit = 3
                ):
        """initialization of ContinuumFit

        Args: 
            continuum_nonezero_ind: pixels to mask at both ends of the order
            base_order_fit: the polynomial order to fit continuum
            max_order_fit: maximum number of the polynomial order for iterate fitting
            nsigma_continuumfit: sigma clipping threshold: tuple (low, high) or list
            maxiter_continuumfit: the maximum number of iteration to fit continuum

        """

        self.continuum_nonzero_ind = continuum_nonzero_ind
        self.base_order_fit = base_order_fit 
        self.max_order_fit = max_order_fit
        self.nsigma_continuumfit = nsigma_continuumfit 
        self.maxiter_continumfit = maxiter_continuumfit

    def fit_continuum(self, x, y, order=6, fitfunc='legendre'):
        """Fit the continuum using sigma clipping

        Args:
            x: the wavelengths
            y: the log-fluxes
            order: the polynomial order to use
            fitfunc: fitting function (default is legendre)

        Returns:
            The value of the continuum at the wavelengths in x

        """
        if fitfunc==None:
            A = np.vander(x - np.nanmean(x), order+1)
        m = np.ones(len(x), dtype=bool)
        for i in range(self.maxiter_continumfit):
            if fitfunc==None:
                w = np.linalg.solve(np.dot(A[m].T, A[m]), np.dot(A[m].T, y[m]))
                mu = np.dot(A, w)
            elif fitfunc=='legendre':
                f = legendre.Legendre.fit(x[m],y[m],order)
                mu = f(x)
            resid = y - mu
            sigma = np.sqrt(np.nanmedian(resid**2))
            m_new = (resid > - self.nsigma_continuumfit[0]*sigma) & (resid < self.nsigma_continuumfit[1]*sigma)
            if m.sum() == m_new.sum():
                m = m_new
                break
            m = m_new
        return mu, m

    def continuum_rsd(self, rsd, npix=2048, ignore_orders=None):
        """fit continuum for rsd

        Args:
            rsd: raw spectrum detector matrix
            npix: number of pixels
            ignore_order: list of orders not to be evaluated the goodness of the fitting

        Returns:
            pandas DataFrame of continuum
        """
        order_fit_itr = self.base_order_fit
        pixels = np.arange(1,npix+1)
        orders = np.arange(1,len(rsd[0])+1)
        if ignore_orders is None:
            ignore_orders = [1,len(rsd[0])-1]

        continuum_nonzero_beginning, continuum_nonzero_ending = self.continuum_nonzero_ind 

        df_continuum = pd.DataFrame([],columns=pixels,index=orders)
        order = orders[0]
        while order != orders[-1]+1:
            flat_ord = rsd[:,order-1]
            cutind = np.zeros(len(flat_ord),dtype=bool)
            cutind[continuum_nonzero_beginning:-continuum_nonzero_ending] = True # mask pixels at both ends of the order
            cutind[np.where(flat_ord==0)[0]] = False # mask zeros
            useind = (~np.isnan(flat_ord))
            continuum_ord, mask = self.fit_continuum(pixels[useind & cutind], flat_ord[useind & cutind], order=order_fit_itr)
            resid = flat_ord[useind & cutind][mask]-continuum_ord[mask]

            if (order not in ignore_orders) and order_fit_itr<self.max_order_fit and 100*np.std(resid/continuum_ord[mask])>1.5:
                # change fitting parameter(s) to avoid bad fitting
                order_fit_itr += 1
                order = orders[0]
            else:
                continuum_ord_interp = np.interp(pixels,pixels[useind & cutind],continuum_ord)
                df_continuum.loc[order,pixels] = continuum_ord_interp
                order += 1

        if order_fit_itr == self.max_order_fit:
            warnings.warn(
                f"Continuum fitting reached the maximum polynomial order ({self.max_order_fit}).\n"
                "          Please check the fitted continuum — it may indicate overfitting.",
                UserWarning
            )
        else:
            print(f"Continuum successfully fitted with polynomial order = {order_fit_itr}.")


        return df_continuum

    def continuum_oneord(self, wdata, flat, order):
        """fit continuum for one order

        Args:
            wdata: the wavelength calibrated target spectrum
            flat: the wavelength calibrated FLAT
            order: order number to fit

        Returns:
            spectrum and continuum of the order
        """

        continuum_nonzero_beginning, continuum_nonzero_ending = self.continuum_nonzero_ind 

        wdata_tmp = wdata[(wdata['order']==order)]
        wav_tmp = wdata_tmp['wav'].values
        flux_tmp = wdata_tmp['flux'].values

        flat_tmp = flat[(flat['order']==order)]
        flatflux_tmp = flat_tmp['flux'].values

        cutind = np.zeros(len(flat_tmp),dtype=bool)
        cutind[continuum_nonzero_beginning:-continuum_nonzero_ending] = True # mask pixels at both ends of the order
        cutind[np.where(flatflux_tmp == 0)[0]] = False # mask zeros
        useind = (~np.isnan(flatflux_tmp)) & (~np.isnan(flux_tmp))

        # CHECK: fit order, clipping sigma
        continuum, mask = self.fit_continuum(wav_tmp[useind & cutind], flatflux_tmp[useind & cutind], order=self.base_order_fit)
        continuum = np.interp(wav_tmp[useind],wav_tmp[useind & cutind],continuum)
        return wav_tmp, flux_tmp, flatflux_tmp, useind, continuum



