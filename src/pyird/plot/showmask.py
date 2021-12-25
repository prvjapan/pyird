import matplotlib.pyplot as plt
import numpy as np
def show_hotpix(obj,im,onoff=False):
    """show hotpixels

    Args:
       obj: sep object
       im: image or mask
       onoff: onoff display mode (use mask as im)

    """
    from matplotlib.patches import Ellipse
    
    # plot background-subtracted image
    fig, ax = plt.subplots()
    m, s = np.mean(im), np.std(im)
    ax.imshow(im, interpolation='nearest', cmap='gray',
                   vmin=m-s, vmax=m+s, origin='lower')
    
    # plot an ellipse for each object
    for i in range(len(obj)):
        e = Ellipse(xy=(obj['x'][i], obj['y'][i]),
                    width=6*obj['a'][i],
                    height=6*obj['b'][i],
                    angle=obj['theta'][i] * 180. / np.pi)
        e.set_facecolor('none')
        e.set_edgecolor('red')
        ax.add_artist(e)
    plt.show()
