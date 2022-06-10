import matplotlib.pyplot as plt
import numpy as np


def show_profile(xprofile_offset_subtracted, yprofile_offset_subtracted,
                 xprofile_offset_subtracted_model, yprofile_offset_subtracted_model):
    """show detector profiles.

    Args:
        xprofile_offset_subtracted
        yprofile_offset_subtracted
        xprofile_offset_subtracted_model
        yprofile_offset_subtracted_model
    """
    nchan = np.shape(xprofile_offset_subtracted)[1]
    fig = plt.figure(figsize=(10, 9))
    ax = fig.add_subplot(311)
    for ichan in range(0, nchan):
        ax.plot(xprofile_offset_subtracted[0, ichan, :], alpha=0.1, color='C0')
        if ichan==0:
            ax.plot(xprofile_offset_subtracted_model[0, :],alpha=1.0, color='C0', label='even pixels')
        else:
            ax.plot(xprofile_offset_subtracted_model[0, :],alpha=1.0, color='C0')

    for ichan in range(0, nchan):
        ax.plot(xprofile_offset_subtracted[1, ichan, :], alpha=0.1, color='C1')
        if ichan==0:
            ax.plot(xprofile_offset_subtracted_model[1, :], alpha=1.0, color='C1', label='odd pixels')
        else:
            ax.plot(xprofile_offset_subtracted_model[1, :], alpha=1.0, color='C1')
            
    ax.set_ylabel('X-profile')
    plt.legend()

    ax2 = fig.add_subplot(312)
    for ichan in range(0, nchan):
        ax2.plot(yprofile_offset_subtracted[0,
                 ichan, :], alpha=0.1, color='C0')
        if ichan==0:
            ax2.plot(yprofile_offset_subtracted_model[0, :], alpha=1.0, color='C0', label='even pixels')
        else:
            ax2.plot(yprofile_offset_subtracted_model[0, :], alpha=1.0, color='C0')

    for ichan in range(0, nchan):
        ax2.plot(yprofile_offset_subtracted[1,
                 ichan, :], alpha=0.1, color='C1')
        if ichan==0:
            ax2.plot(yprofile_offset_subtracted_model[1, :], alpha=1.0, color='C1', label='odd pixels')
        else:
            ax2.plot(yprofile_offset_subtracted_model[1, :], alpha=1.0, color='C1')
            
    ax2.set_ylabel('Y-profile')
    ax3 = fig.add_subplot(313)
    for ichan in range(0, nchan):
        ax3.plot(yprofile_offset_subtracted[0,
                 ichan, :], alpha=0.1, color='C0')
        if ichan==0:
            ax3.plot(yprofile_offset_subtracted_model[0, :], alpha=1.0, color='C0', label='even pixels')
        else:
            ax3.plot(yprofile_offset_subtracted_model[0, :], alpha=1.0, color='C0')
            
    for ichan in range(0, nchan):
        ax3.plot(yprofile_offset_subtracted[1,
                 ichan, :], alpha=0.1, color='C1')
        if ichan==0:
            ax3.plot(yprofile_offset_subtracted_model[1, :], alpha=1.0, color='C1', label='odd pixels')
        else:
            ax3.plot(yprofile_offset_subtracted_model[1, :], alpha=1.0, color='C1')

    ax3.set_ylabel('Y-profile')
    plt.xlim(100, 140)
    plt.ylim(-3, 3)
    plt.xlabel('pixel')
    plt.savefig('profile.png')
    plt.show()


def corrected_detector(im, model_im, corrected_im, vmax=10.0, vmin=-15.0):
    """plot detector images

    Args:
       im: raw image
       model_im: model image
       corrected_im: image after correction
       vmin: color vmin
       vmax: color vmax

    """
    
    fig = plt.figure(figsize=(8, 4))
    ax1 = fig.add_subplot(131)
    cc = ax1.imshow(im, vmin=vmin, vmax=vmax)
    plt.colorbar(cc, shrink=0.55)
    ax1.set_title('Raw image')

    ax2 = fig.add_subplot(132)
    cc = ax2.imshow(model_im, vmin=vmin, vmax=vmax)
    plt.colorbar(cc, shrink=0.55)
    ax2.set_title('Model image')

    ax3 = fig.add_subplot(133)
    cc = ax3.imshow(corrected_im, vmin=-3.0, vmax=4.0)
    plt.colorbar(cc, shrink=0.55)
    ax3.set_title('pattern corrected')
    plt.savefig('pattern_correction.png')
    plt.show()
