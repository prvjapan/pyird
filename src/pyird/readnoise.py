import numpy as np
def RBF(obst,tau):
    if np.shape(np.shape(obst))[0]==1:
        Dt = obst - np.array([obst]).T
    elif np.shape(np.shape(obst))[0]==2:
        Dt = obst
    K=np.exp(-(Dt)**2/2/(tau**2))
    return K

def Matern32(obst,tau):
    if np.shape(np.shape(obst))[0]==1:
        Dt = obst - np.array([obst]).T
    elif np.shape(np.shape(obst))[0]==2:
        Dt = obst
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
    P = invL*(np.dot(Uy.T,np.dot(Dmat,Ux)))
    Kxp=RBFcross(x,xp,xscale)
    Kyp=RBFcross(y,yp,yscale)    
    Dest = (Kxp.T@Ux@P@Uy.T@Kyp)
    return Dest

def RNestimate_OEGP(img,xscale=10,yscale=5,sigma=0.01):
    #odd/even GP method by H.Kawahara
    from numpy import linalg as LA
    import scipy
    import numpy as np
    import time
    #####################
    #GP model for o/e pix and o/e chan
    nchan=32 # number od channel
    noechan=int(nchan/2) #number of odd or even channel
    nRN=64 # number of yaxis
    npix=2048 # number of xaxis (pixel)
    nSRN=int(nRN/2)
    ##fold
    subcube=[]
    for jchan in range(0,nchan):
        arr=img[jchan*nRN:(jchan+1)*nRN,:]
        if np.mod(jchan,2)==0:
            subcube.append(arr[0::2,:])
            subcube.append(arr[1::2,:])
        else:
            subcube.append(arr[-1::-2,:])
            subcube.append(arr[-2::-2,:])
        
    #GP
    subcube_median=np.nanmedian(subcube,axis=(1,2))
    subcubex=subcube-subcube_median[:,np.newaxis,np.newaxis]
    subarray=np.nanmedian(subcubex,axis=0)
    subarray=subarray[:,4:-4]#remove edges

    rectmp=GP2D(subarray,sigma,xscale,yscale)
    rec=np.zeros((nSRN,npix))
    rec[:,4:-4]=rectmp

    ### RECONSTRUCT
    recimg=np.zeros(np.shape(img))
    for jchan in range(0,nchan):
        arr=np.zeros((nRN,npix))
        arr[0::2,:]=rec[:,:]+subcube_median[2*jchan]#+3
        arr[1::2,:]=rec[:,:]+subcube_median[2*jchan+1]
        if np.mod(jchan,2)==0:
            recimg[jchan*nRN:(jchan+1)*nRN,:]=arr
        else:
            recimg[jchan*nRN:(jchan+1)*nRN,:]=arr[::-1,:]
            
    return recimg
#    cimg=imgorig-recimg
