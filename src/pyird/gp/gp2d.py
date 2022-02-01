"""Gaussian Process 2D."""
import numpy as np
from numpy import linalg as LA
import GPy
import matplotlib.pyplot as plt

def GP2Dgpy(X,Z,Nx,Ny):
    """GP 2D for matrix input with mask

    Args:
       X: 2D coordinate (N x 2)
       Z: value (N)
       Nx: x-dimension for a prediction matrix
       Ny: y-dimension for a prediction matrix

    Returns:
       X1,X2,GP prediction (50%)

    """

    xgrid = np.linspace(0, Nx, Nx)
    ygrid = np.linspace(0, Ny, Ny)

    kernel = GPy.kern.Matern32(2, ARD=True)
    model = GPy.models.GPRegression(X, Z, kernel)
    model.optimize(messages=True, max_iters=1e5)

    X1, X2 = np.meshgrid(xgrid, ygrid)
    input_grid = np.array([X1.flatten(), X2.flatten()]).T
    z_pred = model.predict_quantiles(input_grid, quantiles=(50,))[0]
    Zpred=z_pred.reshape(Ny,Nx)

    return X1,X2,Zpred


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
    """GP 2D for different size between input and prediction.

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


if __name__=='__main__':
    np.random.seed(seed=1)

    #data
    Nx = 32
    Ny = 32
    N=Nx*Ny
    xgrid = np.linspace(0, Nx, Nx)
    ygrid = np.linspace(0, Ny, Ny)
    X = (np.array(np.meshgrid(xgrid, ygrid))).reshape(2,Nx*Ny).T
    print(np.shape(X))
    Z = np.sin(X[:, 0:1]/10) * np.sin(X[:, 1:2]/10) + np.random.randn(N, 1)*0.05
       
    # prediction
    X1,X2,Zpred=GP2Dgpy(X,Z,Nx,Ny)
    
    from mpl_toolkits.mplot3d import Axes3D
    fig=plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X1, X2, Zpred, cmap='bwr', linewidth=0,alpha=0.3)
    ax.scatter3D(X[:,0],X[:,1],Z,color="k")
    fig.colorbar(surf)
    plt.show()
    
