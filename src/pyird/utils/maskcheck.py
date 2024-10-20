import warnings


def plotmask(maskfits, obj, vmin=None, vmax=None):
    warn_msg = " Use `plot.showmask.show_mask` instead"
    warnings.warn(warn_msg, FutureWarning)
    from pyird.plot.showmask import show_mask
    show_mask(maskfits, obj, vmin=vmin, vmax=vmax)
