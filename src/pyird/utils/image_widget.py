import astropy.io.fits as pyf
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.cm as cm
import matplotlib.colors as colors

import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import matplotlib
matplotlib.use('tkagg')

class image_1Dand2D(ttk.Frame):
    def __init__(self, master, order, band):
        super().__init__(master)
        self.order = order
        self.band = band

        #self.next_botton()

    def next_botton(self):
        """place a botton on the window"""
        button = tk.Button(self, text="Next / Exit", command=lambda :self.quit_me())
        button.pack()

    def quit_me(self):
        """quit and delete the window"""
        self.quit()
        self.destroy()

    def draw_canvas(self,fig):
        """set the figure on the canvas"""
        ## Show figure on canvas
        canvas = FigureCanvasTkAgg(fig, self)
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        ## Show toolbar
        toolbar=NavigationToolbar2Tk(canvas, self)

        self.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        return canvas

    def title_spec_to_image(self):
        self.winfo_toplevel().title("%s band, Order %d; Press any key anywhere in the spectrum panel!"%(self.band,self.order))

    def title_emission_position(self):
        self.winfo_toplevel().title("%s band, Order %d; Red 'x's are emission like signal, Blue 'x's are bad pixel."%(self.band,self.order))

    def show_spec_to_image(self,rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,wavcal_path=None,hotpix_mask=None,**kwargs):
        """figures of 1d spectrum and 2d detector image

        Args:
            rsd: Raw Spectral Detector matrix
            wav: wavelengths matrix
            mask: mask of trace lines
            pixcoord, iys_plot, iye_plot: pixel coordinate of x and y direction
            rotim: rotated image (intend to align trace orientation)
            wavcal_path: path to master file of the wavelength calibration
            hotpix_mask: hotpixel mask
        """

        if not any(kwargs):
            kwargs = {'vmin':-10,'vmax':50,'scale':'linear'}

        fig=plt.figure()#figsize=(9,6))
        ax1=fig.add_subplot(211)

        ## 1D spectrum
        rsd_ord = np.nan_to_num(rsd[:,self.order-1],0)
        if wavcal_path is None: # x axis will be 'pixels'
            ax1.plot(rsd_ord)
            text=ax1.text(0,0, "", va="bottom", ha="left")
            ax1.set(xlabel='pixel')
        else:
            wav_ord = wav[:,self.order-1]
            ax1.plot(wav_ord,rsd_ord)
            text=ax1.text(min(wav_ord),0, "", va="bottom", ha="left")
            ax1.set(xlabel='wavelength [nm]',xlim=(min(wav_ord),max(wav_ord)))

        ## 2D detector image
        ax2=fig.add_subplot(212)
        if kwargs['scale']=='log':
            rotim_cp = rotim.copy()
            rotim_cp[rotim_cp<0] = 1
            im_plot = np.log(rotim_cp.T)
        elif kwargs['scale']=='linear':
            im_plot = rotim.T
        image = ax2.imshow(im_plot,cmap="OrRd",vmin=kwargs['vmin'],vmax=kwargs['vmax'],origin="lower")
        if not hotpix_mask is None:
            hpmask_plot = hotpix_mask.astype(int)*kwargs['vmax']
            #hpmask_plot[~hotpix_mask] = kwargs['vmin']
            hpmask_plot = hpmask_plot[::-1,::-1].T#.T[::-1,::-1]
            #show hotpix as imshow
            cmap = cm.cool_r
            cmap_data = cmap(np.arange(cmap.N))
            cmap_data[0, 3] = 0 # 0 のときのα値を0(透明)にする
            customized_cool = colors.ListedColormap(cmap_data)
            ax2.imshow(hpmask_plot,cmap=customized_cool,origin="lower",alpha=0.9)

        """
            #show hotpix as contour #CAUTION!! TAKES LONG TIME TO PLOT (~a few tens of second)
            fac=10
            condition = np.kron(hpmask_plot > 0, np.ones((fac, fac)))
            extent = (-0.5, hpmask_plot.shape[1]-0.5, -0.5, hpmask_plot.shape[0]-0.5)
            ax2.contour(condition, levels=[0.5], extent=extent,cmap='cool',linewidths=[1])
        """

        ## show traced apertures
        for i in range(len(pixcoord)):
            ax2.fill_between(pixcoord[i],iys_plot[i],iye_plot[i],color='tab:blue',alpha=0.5)
        fig.colorbar(image)

        ## plot 'x' when a key is pressed
        def onkey(event):
            xi = event.xdata
            cmap = plt.cm.get_cmap("Set1",9)
            color = cmap(int(xi*1e2)%9)
            rsd_ord = np.nan_to_num(rsd[:,self.order-1],0)
            if wavcal_path is not None:
                wav_ord = wav[:,self.order-1]
                xi = np.where(np.abs(wav_ord-float(xi))==min(np.abs(wav_ord-float(xi))))[0][0] #'searchsorted' doesn't work
                tx = 'key=%s, xdata=%.3f, order=%d, wavcal=%s, pix=%.1f' % (event.key, event.xdata,self.order,bool(wavcal_path),int(xi))
                ax1.scatter(wav_ord[xi],rsd_ord[xi],color=color,marker='x',s=5,lw=12)
            else:
                tx = 'key=%s, xdata=%.3f, order=%d, wavcal=%s' % (event.key, event.xdata,self.order,bool(wavcal_path))
                ax1.scatter(int(xi),rsd_ord[xi],color=color,marker='x',s=5,lw=12)
            ax2.scatter(int(xi),mask[self.order-1][int(xi)],color=color,marker='x',s=5,lw=12)
            text.set_text(tx)
            canvas.draw()

        cid = fig.canvas.mpl_connect('key_press_event', onkey)

        ax1.set(ylabel='flux')#,ylim=(-0.01,100))

        canvas = self.draw_canvas(fig)
        self.title_spec_to_image()

    def show_emission_position(self,stream2D,rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,wavcal_path=None,hotpix_mask=None,fit=True,**kwargs):
        """detect emissions on the spectrum and detector image of an arbitral aperture.

        Args:
            Stream2D: pyird.utils.stream2D (need to have trace information)
            rsd: Raw Spectral Detector matrix
            wav: wavelengths matrix
            mask: mask of trace lines
            pixcoord, iys_plot, iye_plot: pixel coordinate of x and y direction
            rotim: rotated image (intend to align trace orientation)
            wavcal_path: path to master file of the wavelength calibration
            hotpix_mask: hotpixel mask
        """
        xmin, xmax = stream2D.trace.xmin, stream2D.trace.xmax
        npix = 2048
        if not any(kwargs):
            kwargs = {'vmin':-10,'vmax':50,'scale':'linear'}

        ## make images of each order (npix x aperture width)
        j=self.order-1
        ap_order=[]
        ap_order.extend(np.zeros((xmin[j],iye_plot[j][0]-iys_plot[j][0])))
        for i in range(len(pixcoord[j])):
            ap_order.append(rotim[pixcoord[j][i],iys_plot[j][i]:iye_plot[j][i]])
        ap_order.extend(np.zeros((npix-1-xmax[j],iye_plot[j][0]-iys_plot[j][0])))
        ap_order=np.array(ap_order).T

        if not hotpix_mask is None:
            hpmask = hotpix_mask.astype(int)
            hpmask[hotpix_mask] = kwargs['vmax']
            hpmask = hpmask[::-1,::-1]
            hpmask_order=[]
            hpmask_order.extend(np.zeros((xmin[j],iye_plot[j][0]-iys_plot[j][0])))
            for i in range(len(pixcoord[j])):
                hpmask_order.append(hpmask[pixcoord[j][i],iys_plot[j][i]:iye_plot[j][i]])
            hpmask_order.extend(np.zeros((npix-1-xmax[j],iye_plot[j][0]-iys_plot[j][0])))
            hpmask_order = np.array(hpmask_order).T

        ## detection of emission on 2d image
        from scipy.ndimage.filters import maximum_filter
        def detect_peaks(image, filter_size=3, order=0.5):
            local_max = maximum_filter(image, footprint=np.ones((filter_size, filter_size)), mode='constant')
            detected_peaks = np.ma.array(image, mask=~(image == local_max))

            ## exclude peaks lower than 'order'
            temp = np.ma.array(detected_peaks, mask=~(detected_peaks >= order))#detected_peaks.max() * order))
            peaks_index = np.where((temp.mask != True))
            return peaks_index

        maxid = detect_peaks(ap_order, filter_size=5, order=max(3,np.median(ap_order)*2)) ##CHECK!!

        ## plot image (npix x aperture width)
        fig=plt.figure()#figsize=(9,6))
        ax1=fig.add_subplot(211)
        #ax1.imshow(ap_order,aspect=10,cmap="OrRd",vmin=-10,vmax=50,origin="lower")
        if kwargs['scale']=='log':
            ap_cp = ap_order.copy()
            ap_cp[ap_cp<=0] = 1
            ap_order_plot = np.log(ap_cp)
        elif kwargs['scale']=='linear':
            ap_order_plot = ap_order
        image = ax1.imshow(ap_order_plot,aspect=10,cmap="OrRd",vmin=kwargs['vmin'],vmax=np.max(ap_order_plot),origin="lower")
        if not hotpix_mask is None:
            fac=10
            condition = np.kron(hpmask_order > 0, np.ones((fac, fac)))
            extent = (-0.5, hpmask_order.shape[1]-0.5, -0.5, hpmask_order.shape[0]-0.5)
            ax1.contour(condition, levels=[0.5], extent=extent,cmap='cool',linewidths=[1])

        fig.colorbar(image)

        if fit:
            from astropy.modeling import models, fitting
            amp=1000
            sigma=1
            fit_w = fitting.LevMarLSQFitter()
            skylike_ord = []
            for i in range(len(maxid[0])):
                #warnings.resetwarnings()
                w = models.Gaussian2D(amp, maxid[1][i], maxid[0][i], sigma, sigma)
                yi, xi = np.indices(ap_order.shape)
                g = fit_w(w, xi, yi, ap_order)
                ax1.scatter(maxid[1][i], maxid[0][i], color='skyblue',marker='x',alpha=0.5)
                if not fit_w.fit_info['message'].startswith('Both'):
                    #print(fit_w.fit_info['message'])
                    continue
                model_data = g(xi, yi)

                if (g.y_mean.value<1) or (len(ap_order)-1<g.y_mean.value) or ((g.x_stddev.value<1) and (g.y_stddev.value<1)) or (np.abs(maxid[1][i]-g.x_mean.value)>1.5):
                    ax1.scatter(maxid[1][i], maxid[0][i], color='blue',marker='x')
                    #circle = patches.Circle((g.x_mean.value, g.y_mean.value),
                    #         g.x_stddev.value, facecolor ='none',
                    #        edgecolor = 'skyblue', linewidth = 1)
                    #ax.add_patch(circle)
                else:
                    ax1.scatter(maxid[1][i], maxid[0][i], color='red',marker='x')
                    circle = patches.Circle((g.x_mean.value, g.y_mean.value),
                                     g.x_stddev.value, facecolor ='none',
                                    edgecolor = 'pink', linewidth = 1)
                    ax1.add_patch(circle)
                    skylike_ord.append(g.x_mean.value)

        ## plot spectrum of an order
        ax2=fig.add_subplot(212)
        ax3=ax2.twiny()
        if wavcal_path is not None:
            ax2.plot(wav[:,j],np.nan_to_num(rsd[:,j],0),alpha=0.5)
            #for i in range(len(skylike_ord)):
            #    ax2.vlines(wav[int(skylike_ord[i]),j],0,10000,color='pink',alpha=0.5)
            ax2.set(xlabel='wavelength [nm]',xlim=(min(wav[:,j]),max(wav[:,j])))

            #pix2wav = lambda x: wav[x,j]
            wav2pix = lambda x: np.arange(0,npix+1,1)[np.searchsorted(wav[:,j],x)]
            xvals = ax2.get_xticks()
            ax3.plot(wav[:,j],np.nan_to_num(rsd[:,j],0),alpha=0.1)
            ax3.set_xticks(xvals,labels=wav2pix(xvals))
            ax3.set(xlabel='pixel',xlim=(min(wav[:,j]),max(wav[:,j])))
        else:
            ax2.plot(np.nan_to_num(rsd[:,j],0),alpha=0.5)
            ax2.set(xlabel='pixel')
        ax2.set(ylabel='flux')

        canvas = self.draw_canvas(fig)
        self.title_emission_position()
