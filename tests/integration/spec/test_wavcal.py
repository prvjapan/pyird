import pytest
from pyird.spec.wavcal import fit_polynomial, fit_wav_solution

def test_wavcal():
    import numpy as np
    npix = 2048
    j, l = 0, 21
    pixels = np.arange(1,npix+1,1)
    orders = np.arange(1,l+1-j,1)
    X,Y = np.meshgrid(pixels,orders)
    Ni, Nx = 5, 4

    #example of coefficients
    p_in = [[5.58730004e+00, -1.71104223e-03,  7.50009100e-07, -9.94465893e-11], \
             [1.39281235e+03,  1.39706356e-02, -1.75209173e-06,  7.01238665e-11], \
             [9.00233492e+00, -8.57053225e-05,  7.97588392e-08, -1.46187772e-11], \
             [4.13355786e-02,  8.28862001e-06, -4.62958089e-09,  8.39980739e-13],\
             [5.08278747e-04, -1.25877266e-07,  7.84194051e-11, -1.51053784e-14]]

    w = np.ones(X.shape)

    data_in = fit_polynomial((X,Y),Ni,Nx,params=p_in).reshape(npix,l-j)
    p_out = fit_wav_solution((X,Y),data_in,w,Ni,Nx)
    p_relerr = (p_out/p_in) - 1
    #data_out = fitfunc((X,Y),Ni,Nx,p_out).reshape(npix,l-j)
    assert np.max(p_relerr) < 1e-8

if __name__ == '__main__':
    test_wavcal()
