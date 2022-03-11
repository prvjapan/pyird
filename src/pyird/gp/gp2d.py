"""Gaussian Process 2D.

- by defalt, pyird uses gpkron package for a 2D GP. However, this module uses GPy instead.


"""


import numpy as np
from numpy import linalg as LA
import GPy
import matplotlib.pyplot as plt


def GP2Dgpy(X, Z, Nx, Ny):
    """GP 2D for matrix input with mask.

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
    Zpred = z_pred.reshape(Ny, Nx)

    return X1, X2, Zpred


if __name__ == '__main__':
    np.random.seed(seed=1)

    # data
    Nx = 32
    Ny = 32
    N = Nx*Ny
    xgrid = np.linspace(0, Nx, Nx)
    ygrid = np.linspace(0, Ny, Ny)
    X = (np.array(np.meshgrid(xgrid, ygrid))).reshape(2, Nx*Ny).T
    print(np.shape(X))
    Z = np.sin(X[:, 0:1]/10) * np.sin(X[:, 1:2]/10) + \
        np.random.randn(N, 1)*0.05

    # prediction
    X1, X2, Zpred = GP2Dgpy(X, Z, Nx, Ny)

    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X1, X2, Zpred, cmap='bwr', linewidth=0, alpha=0.3)
    ax.scatter3D(X[:, 0], X[:, 1], Z, color='k')
    fig.colorbar(surf)
    plt.show()
