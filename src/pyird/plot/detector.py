import matplotlib.pyplot as plt
import numpy as np

def show_profile(xprofile_offset_subtracted,yprofile_offset_subtracted,\
                 xprofile_offset_subtracted_model,yprofile_offset_subtracted_model):
    """show detector profiles

    Args:
        xprofile_offset_subtracted
        yprofile_offset_subtracted
        xprofile_offset_subtracted_model
        yprofile_offset_subtracted_model
    
    """
    nchan=np.shape(xprofile_offset_subtracted)[1]
    fig=plt.figure(figsize=(10,9))
    ax=fig.add_subplot(311)
    for ichan in range(0,nchan):
        ax.plot(xprofile_offset_subtracted[0,ichan,:],alpha=0.1, color="C0")
        ax.plot(xprofile_offset_subtracted_model[0,:],alpha=1.0, color="C0",label="even pixels")
    for ichan in range(0,nchan):
        ax.plot(xprofile_offset_subtracted[1,ichan,:],alpha=0.1, color="C1")
        ax.plot(xprofile_offset_subtracted_model[1,:],alpha=1.0, color="C1",label="odd pixels")
    ax.set_ylabel("X-profile")
    plt.legend()

    ax2=fig.add_subplot(312)
    for ichan in range(0,nchan):
        ax2.plot(yprofile_offset_subtracted[0,ichan,:],alpha=0.1, color="C0")
        ax2.plot(yprofile_offset_subtracted_model[0,:],alpha=1.0, color="C0",label="even pixels")
    for ichan in range(0,nchan):
        ax2.plot(yprofile_offset_subtracted[1,ichan,:],alpha=0.1, color="C1")
        ax2.plot(yprofile_offset_subtracted_model[1,:],alpha=1.0, color="C1",label="odd pixels")
    ax2.set_ylabel("Y-profile")
    ax3=fig.add_subplot(313)
    for ichan in range(0,nchan):
        ax3.plot(yprofile_offset_subtracted[0,ichan,:],alpha=0.1, color="C0")
        ax3.plot(yprofile_offset_subtracted_model[0,:],alpha=1.0, color="C0",label="even pixels")
    for ichan in range(0,nchan):
        ax3.plot(yprofile_offset_subtracted[1,ichan,:],alpha=0.1, color="C1")
        ax3.plot(yprofile_offset_subtracted_model[1,:],alpha=1.0, color="C1",label="odd pixels")
    ax3.set_ylabel("Y-profile")
    plt.xlim(100,140)
    plt.ylim(-3,3)
    plt.xlabel("pixel")
    plt.savefig("profile.png")
    plt.show()
