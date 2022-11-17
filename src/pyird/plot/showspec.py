import matplotlib.pyplot as plt
from pyird.image.hotpix import apply_hotpixel_mask

def show_wavcal_spectrum(df_wspec,**kwargs):
    """show the wavelength calibrated spectrum.

    Args:
        df_wspec: pandas.DataFrame that contains wavelengths and flux
        kwargs: keys for plot

    """
    wav = df_wspec['wav'].values
    flux_col = df_wspec.columns[[('flux' in x) for x in df_wspec.columns]]
    flux = df_wspec[flux_col].values

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(wav,flux,**kwargs)
    plt.show()
"""
def show_spec_to_image(im,hotpix_mask=None,wavcal_path=None):

    rawspec, pixcoord, rotim, tl, iys_plot, iye_plot = flatten(
        im, trace_legendre, y0, xmin, xmax, coeff, self.inst)

    mask_shape = (2048,2048)
    mask = np.zeros(mask_shape)
    for i in range(len(y0)):
        mask[i,xmin[i]:xmax[i]+1] = tl[i]

    rsd = multiorder_to_rsd(rawspec, pixcoord)

    if not hotpix_mask is None:
        rsd = apply_hotpixel_mask(hotpix_mask, rsd, y0, xmin, xmax, coeff)

    if not wavcal_path is None:
        #master_path = thar.anadir/('thar_%s_%s.fits'%(thar.band,thar.trace.mmf)) ##CHECK!!
        hdu = pyf.open(wavcal_path)
        wav = hdu[0].data
        hdu.close()
    else:
        wav = []

    fig=plt.figure(figsize=(9,6))
    def plot_rsd(order):
        fig.clear()
        ax1=fig.add_subplot(211)
        if wavcal==False:
            #ax.plot(rsd[:,order-1])
            ax1.plot(np.nan_to_num(rsd[:,order-1],0))
            text=ax1.text(0,0, "", va="bottom", ha="left")
            ax1.set(xlabel='pixel')
        else:
            wav_tmp = wav[:,order-1]
            ax1.plot(wav_tmp,np.nan_to_num(rsd[:,order-1],0))
            text=ax1.text(min(wav[:,order-1]),0, "", va="bottom", ha="left")
            ax1.set(xlabel='wavelength [nm]',xlim=(min(wav[:,order-1]),max(wav[:,order-1])))

        ax2=fig.add_subplot(212)
        if kwargs['scale'].value=='log':
            rotim_cp = rotim.copy()
            rotim_cp[rotim_cp<0] = 1
            im_plot = np.log(rotim_cp.T)
        elif kwargs['scale'].value=='linear':
            im_plot = rotim.T
        image = ax2.imshow(im_plot,cmap="OrRd",vmin=kwargs['vmin'].value,vmax=kwargs['vmax'].value,origin="lower")
        if not hotpix_mask is None:
            hpmask_plot = hotpix_mask.astype(int)*kwargs['vmax'].value
            #hpmask_plot[~hotpix_mask] = kwargs['vmin'].value
            hpmask_plot = hpmask_plot[::-1,::-1].T#.T[::-1,::-1]
            #show hotpix as imshow
            cmap = cm.cool_r
            cmap_data = cmap(np.arange(cmap.N))
            cmap_data[0, 3] = 0 # 0 のときのα値を0(透明)にする
            customized_cool = colors.ListedColormap(cmap_data)
            ax2.imshow(hpmask_plot,cmap=customized_cool,origin="lower",alpha=0.9)

        for i in range(len(y0)):
            ax2.fill_between(pixcoord[i],iys_plot[i],iye_plot[i],color='tab:blue',alpha=0.5)
        fig.colorbar(image)

        def onkey(event):
            xi = event.xdata
            if wavcal==True:
                tx = 'key=%s, xdata=%f, order=%d, wavcal=%s, pix=%f' % (event.key, event.xdata,order,wavcal,np.searchsorted(wav[:,order-1],event.xdata))
                wav_tmp = wav[:,order-1]
                xi = np.where(np.abs(wav_tmp-xi)==min(np.abs(wav_tmp-xi)))[0][0] #searchsortedだとうまくいかない
            else:
                tx = 'key=%s, xdata=%f, order=%d, wavcal=%s' % (event.key, event.xdata,order,wavcal)
            ax2.plot(int(xi),mask[order-1][int(xi)],'rx')
            text.set_text(tx)

        cid = fig.canvas.mpl_connect('key_press_event', onkey)

        ax1.set(ylabel='flux')#,ylim=(-0.01,100))

    #interact(plot_rsd,order=range(1,rsd.shape[1]+1))
    plot_rsd(order=10)

    #fig.canvas.mpl_disconnect(cid)
    return rotim,pixcoord,iys_plot,iye_plot,xmin,xmax,wav,rsd
"""
