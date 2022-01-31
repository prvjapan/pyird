"""Gaussian Process 2D


"""

from numpy import linalg as LA
import numpy as np


def RBF(obst, tau):
    """RBF kernel
    Args:
        obst: variable
        tau: scale
    Returns
        K
    """
    Dt = obst - np.array([obst]).T
    K = np.exp(-(Dt)**2/2/(tau**2))
    return K


def RBFcross(obst, pret, tau):
    Dt = obst - np.array([pret]).T
    K = np.exp(-(Dt)**2/2/(tau**2))
    return K


def Matern32(obst, tau):
    Dt = obst - np.array([obst]).T
    fac = np.sqrt(3.0)*np.abs(Dt)/tau
    K = (1.0+fac)*np.exp(-fac)
    return K


def Matern32cross(obst, pret, tau):
    Dt = obst - np.array([pret]).T
    fac = np.sqrt(3.0)*np.abs(Dt)/tau
    K = (1.0+fac)*np.exp(-fac)
    return K


def GP2D(Dmat, sigma, xscale, yscale):
    Nx, Ny = np.shape(Dmat)
    x = np.array(list(range(0, Nx)))
    y = np.array(list(range(0, Ny)))
    Kx = Matern32(x, xscale)
    Ky = Matern32(y, yscale)
    kapx, Ux = LA.eigh(Kx)
    np.diag(kapx)
    kapy, Uy = LA.eigh(Ky)
    np.diag(kapy)
    invL = 1.0/(np.outer(kapx, kapy)+sigma**2)
    P = invL*(np.dot(Ux.T, np.dot(Dmat, Uy)))
    #    P = invL*(Ux.T@Dmat@Uy)
    Dest = Dmat - (sigma**2)*(np.dot(Ux, np.dot(P, Uy.T)))
    #    Dest = Dmat - (sigma**2)*(Ux@P@Uy.T)
    return Dest


def GP2Dcross(Dmat, Dpre, sigma, xscale, yscale):
    """GP 2D for different size between input and prediction

    Note:
       It should be Dmat <= Dpre. See #12 for equations.

    Args:
        Dmat: input 2D matrix
        Dpre: prediction matrix form (empty)
        sigma: observational Gaussian noise std
        xscale: GP correlated length (hyperparameter) for X
        yscale: GP correlated length (hyperparameter) for Y

    """
    rat = np.array(np.shape(Dpre))/np.array(np.shape(Dmat))
    Nx, Ny = np.shape(Dmat)

    x = (np.array(list(range(0, Nx)))+0.5)*rat[0]
    y = (np.array(list(range(0, Ny)))+0.5)*rat[1]

    Nxp, Nyp = np.shape(Dpre)
    xp = np.array(list(range(0, Nxp)))
    yp = np.array(list(range(0, Nyp)))

    Kx = RBF(x, xscale)
    Ky = RBF(y, yscale)
    kapx, Ux = LA.eigh(Kx)
    np.diag(kapx)
    kapy, Uy = LA.eigh(Ky)
    np.diag(kapy)
    invL = 1.0/(np.outer(kapx, kapy)+sigma**2)
    P = invL*(np.dot(Ux.T, np.dot(Dmat, Uy)))
    Kxp = RBFcross(x, xp, xscale)
    Kyp = RBFcross(y, yp, yscale)
    Dest = (Kxp@Ux@P@Uy.T@Kyp.T)
    return Dest
