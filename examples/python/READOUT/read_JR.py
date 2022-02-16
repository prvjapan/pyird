import numpy as np
import matplotlib.pyplot as plt
import tqdm
model_im=np.load("readout_pattern.npz")["arr_0"]
l,m,n=np.shape(model_im)
for i in tqdm.tqdm(range(0,l)):
    c=plt.imshow(model_im[i,:,:],vmin=-20.0,vmax=20.0)
    plt.colorbar(c,shrink=0.7)
    plt.savefig("fig/"+str(i)+".png")
    plt.clf()
