from numpy import linalg as LA
import scipy
import numpy as np

def RBF(obst,tau):
    Dt = obst - np.array([obst]).T
    K=np.exp(-(Dt)**2/2/(tau**2))
    return K

def RBFcross(obst,pret,tau):
    Dt = obst - np.array([pret]).T
    K=np.exp(-(Dt)**2/2/(tau**2))
    return K

def Matern32(obst,tau):
    Dt = obst - np.array([obst]).T
    fac=np.sqrt(3.0)*np.abs(Dt)/tau
    K=(1.0+fac)*np.exp(-fac)
    return K

def Matern32cross(obst,pret,tau):
    Dt = pret - np.array([obst]).T
    fac=np.sqrt(3.0)*np.abs(Dt)/tau
    K=(1.0+fac)*np.exp(-fac)
    return K

def GP2D(Dmat,sigma,xscale,yscale):
    Nx,Ny=np.shape(Dmat)
    x=np.array(list(range(0,Nx)))
    y=np.array(list(range(0,Ny)))
    Kx=Matern32(x,xscale)
    Ky=Matern32(y,yscale)
    kapx,Ux=LA.eigh(Kx)
    Lx=np.diag(kapx)
    kapy,Uy=LA.eigh(Ky)
    Ly=np.diag(kapy)
    invL=1.0/(np.outer(kapx,kapy)+sigma**2)
    P = invL*(np.dot(Ux.T,np.dot(Dmat,Uy)))
    #    P = invL*(Ux.T@Dmat@Uy)
    Dest = Dmat - (sigma**2)*(np.dot(Ux,np.dot(P,Uy.T)))
    #    Dest = Dmat - (sigma**2)*(Ux@P@Uy.T)
    return Dest


def GP2Dcross(Dmat,Dpre,sigma,xscale,yscale):
    rat=np.array(np.shape(Dpre))/np.array(np.shape(Dmat))
    Nx,Ny=np.shape(Dmat)
    
    x=(np.array(list(range(0,Nx)))+0.5)*rat[0]
    y=(np.array(list(range(0,Ny)))+0.5)*rat[1]

    Nxp,Nyp=np.shape(Dpre)
    xp=np.array(list(range(0,Nxp)))
    yp=np.array(list(range(0,Nyp)))

    Kx=RBF(x,xscale)
    Ky=RBF(y,yscale)
    kapx,Ux=LA.eigh(Kx)
    Lx=np.diag(kapx)
    kapy,Uy=LA.eigh(Ky)
    Ly=np.diag(kapy)
    invL=1.0/(np.outer(kapx,kapy)+sigma**2)
    P = invL*(np.dot(Ux.T,np.dot(Dmat,Uy)))
    Kxp=RBFcross(x,xp,xscale)
    Kyp=RBFcross(y,yp,yscale)    
    Dest = (Kxp@Ux@P@Uy.T@Kyp.T)
    return Dest
