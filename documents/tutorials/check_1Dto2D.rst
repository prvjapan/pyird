Tools for Mapping 1D Spectra to 2D Detector Images
==================================================

Sometimes we want to check the detector image corresponds to a point of
the extracted spectrum.

Preparations
~~~~~~~~~~~~

The first part is the same as the usual reduction process before assign
wavelength to the extracted spectrum. See the tutorial of IRD_stream for
more detail.

.. code:: ipython3

    from pyird.utils import irdstream
    import pathlib
    from pyird.image.bias import bias_subtract_image
    from pyird.image.hotpix import identify_hotpix_sigclip
    import astropy.io.fits as pyf
    
    #--------SETTINGS--------#
    basedir = pathlib.Path('~/pyird/data/20210317/').expanduser()
    
    band = 'h' #'h' or 'y'
    mmf = 'mmf2' #'mmf1' (comb fiber) or 'mmf2' (star fiber)
    readout_noise_mode = "real" #'real' or 'default'
    
    datadir_flat = basedir/'flat/'
    datadir_dark = basedir/'dark/'
    datadir_thar = basedir/'thar'
    datadir_target = basedir/'target/'
    anadir = basedir/'reduc/'
    
    fitsid_flat_comb = list(range(41704,41804,2)) 
    fitsid_flat_star = (range(41804,41904,2)) 
    fitsid_dark = [41504]
    fitsid_thar = list(range(14632,14732))
    fitsid_target = [41510]
    #-------------------------#
    
    #--------FOR CALIBRATION--------#
    # aperture extraction
    flat=irdstream.Stream2D("flat",datadir_flat,anadir)
    flat.fitsid=fitsid_flat_comb ##FLAT_COMB
    
    ################################
    ### SELECT H band or YJ band ###
    ################################
    flat.band=band
    if flat.band=='h':
        flat.fitsid_increment() # when you use H-band
        trace_mmf=flat.aptrace(cutrow = 800,nap=42) #TraceAperture instance
    elif flat.band=='y':
        trace_mmf=flat.aptrace(cutrow = 1000,nap=102) #TraceAperture instance
    
    # hotpixel mask
    dark = irdstream.Stream2D('dark', datadir_dark, anadir,fitsid=fitsid_dark)
    if flat.band=='h':
        dark.fitsid_increment() # when you use H-band
    for data in dark.rawpath:
        im = pyf.open(str(data))[0].data
    im_subbias = bias_subtract_image(im)
    hotpix_mask = identify_hotpix_sigclip(im_subbias)
    
    if mmf=="mmf2":
        trace_mmf.choose_mmf2_aperture() #mmf2 (star fiber)
    elif mmf=="mmf1":
        trace_mmf.choose_mmf1_aperture() #mmf1 (comb fiber)
    
    # load ThAr raw image
    if band=='h':
        rawtag='IRDAD000'
    elif band=='y':
        rawtag='IRDBD000'
    
    #wavelength calibration
    thar=irdstream.Stream2D("thar",datadir_thar,anadir,rawtag=rawtag,fitsid=fitsid_thar)
    thar.trace = trace_mmf
    thar.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)
    thar.calibrate_wavelength()
    
    ### TARGET ###
    # Load data
    target = irdstream.Stream2D('targets', datadir_target, anadir, fitsid=fitsid_target)
    if flat.band=='h':
        target.fitsid_increment() # when you use H-band
    target.info = True  # show detailed info
    target.trace = trace_mmf
    # clean pattern
    target.clean_pattern(extin='', extout='_cp', hotpix_mask=hotpix_mask)



.. parsed-literal::

    No fitsid yet.
    median combine:  


.. parsed-literal::

    100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 50/50 [00:00<00:00, 417.16it/s]


.. parsed-literal::

    default nap value
    cross-section: row  1170



.. image:: check_1Dto2D_files/check_1Dto2D_1_3.png


.. parsed-literal::

    100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 42/42 [00:12<00:00,  3.47it/s]



.. image:: check_1Dto2D_files/check_1Dto2D_1_5.png


.. parsed-literal::

    fitsid: [41504]
    hotpix mask = 0.58 percent
    fitsid: [14632, 14633, 14634, 14635, 14636, 14637, 14638, 14639, 14640, 14641, 14642, 14643, 14644, 14645, 14646, 14647, 14648, 14649, 14650, 14651, 14652, 14653, 14654, 14655, 14656, 14657, 14658, 14659, 14660, 14661, 14662, 14663, 14664, 14665, 14666, 14667, 14668, 14669, 14670, 14671, 14672, 14673, 14674, 14675, 14676, 14677, 14678, 14679, 14680, 14681, 14682, 14683, 14684, 14685, 14686, 14687, 14688, 14689, 14690, 14691, 14692, 14693, 14694, 14695, 14696, 14697, 14698, 14699, 14700, 14701, 14702, 14703, 14704, 14705, 14706, 14707, 14708, 14709, 14710, 14711, 14712, 14713, 14714, 14715, 14716, 14717, 14718, 14719, 14720, 14721, 14722, 14723, 14724, 14725, 14726, 14727, 14728, 14729, 14730, 14731]
    clean_pattern: output extension=_cp


.. parsed-literal::

    100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 21/21 [00:00<00:00, 90.20it/s]


.. parsed-literal::

    Ignore IRDAD00014632.fits -> IRDAD00014632_cp.fits
    Ignore IRDAD00014633.fits -> IRDAD00014633_cp.fits
    Ignore IRDAD00014634.fits -> IRDAD00014634_cp.fits
    Ignore IRDAD00014635.fits -> IRDAD00014635_cp.fits
    Ignore IRDAD00014636.fits -> IRDAD00014636_cp.fits
    Ignore IRDAD00014637.fits -> IRDAD00014637_cp.fits
    Ignore IRDAD00014638.fits -> IRDAD00014638_cp.fits
    Ignore IRDAD00014639.fits -> IRDAD00014639_cp.fits
    Ignore IRDAD00014640.fits -> IRDAD00014640_cp.fits
    Ignore IRDAD00014641.fits -> IRDAD00014641_cp.fits
    Ignore IRDAD00014642.fits -> IRDAD00014642_cp.fits
    Ignore IRDAD00014643.fits -> IRDAD00014643_cp.fits
    Ignore IRDAD00014644.fits -> IRDAD00014644_cp.fits
    Ignore IRDAD00014645.fits -> IRDAD00014645_cp.fits
    Ignore IRDAD00014646.fits -> IRDAD00014646_cp.fits
    Ignore IRDAD00014647.fits -> IRDAD00014647_cp.fits
    Ignore IRDAD00014648.fits -> IRDAD00014648_cp.fits
    Ignore IRDAD00014649.fits -> IRDAD00014649_cp.fits
    Ignore IRDAD00014650.fits -> IRDAD00014650_cp.fits
    Ignore IRDAD00014651.fits -> IRDAD00014651_cp.fits
    Ignore IRDAD00014652.fits -> IRDAD00014652_cp.fits
    Ignore IRDAD00014653.fits -> IRDAD00014653_cp.fits
    Ignore IRDAD00014654.fits -> IRDAD00014654_cp.fits
    Ignore IRDAD00014655.fits -> IRDAD00014655_cp.fits
    Ignore IRDAD00014656.fits -> IRDAD00014656_cp.fits
    Ignore IRDAD00014657.fits -> IRDAD00014657_cp.fits
    Ignore IRDAD00014658.fits -> IRDAD00014658_cp.fits
    Ignore IRDAD00014659.fits -> IRDAD00014659_cp.fits
    Ignore IRDAD00014660.fits -> IRDAD00014660_cp.fits
    Ignore IRDAD00014661.fits -> IRDAD00014661_cp.fits
    Ignore IRDAD00014662.fits -> IRDAD00014662_cp.fits
    Ignore IRDAD00014663.fits -> IRDAD00014663_cp.fits
    Ignore IRDAD00014664.fits -> IRDAD00014664_cp.fits
    Ignore IRDAD00014665.fits -> IRDAD00014665_cp.fits
    Ignore IRDAD00014666.fits -> IRDAD00014666_cp.fits
    Ignore IRDAD00014667.fits -> IRDAD00014667_cp.fits
    Ignore IRDAD00014668.fits -> IRDAD00014668_cp.fits
    Ignore IRDAD00014669.fits -> IRDAD00014669_cp.fits
    Ignore IRDAD00014670.fits -> IRDAD00014670_cp.fits
    Ignore IRDAD00014671.fits -> IRDAD00014671_cp.fits
    Ignore IRDAD00014672.fits -> IRDAD00014672_cp.fits
    Ignore IRDAD00014673.fits -> IRDAD00014673_cp.fits
    Ignore IRDAD00014674.fits -> IRDAD00014674_cp.fits
    Ignore IRDAD00014675.fits -> IRDAD00014675_cp.fits
    Ignore IRDAD00014676.fits -> IRDAD00014676_cp.fits
    Ignore IRDAD00014677.fits -> IRDAD00014677_cp.fits
    Ignore IRDAD00014678.fits -> IRDAD00014678_cp.fits
    Ignore IRDAD00014679.fits -> IRDAD00014679_cp.fits
    Ignore IRDAD00014680.fits -> IRDAD00014680_cp.fits
    Ignore IRDAD00014681.fits -> IRDAD00014681_cp.fits
    Ignore IRDAD00014682.fits -> IRDAD00014682_cp.fits
    Ignore IRDAD00014683.fits -> IRDAD00014683_cp.fits
    Ignore IRDAD00014684.fits -> IRDAD00014684_cp.fits
    Ignore IRDAD00014685.fits -> IRDAD00014685_cp.fits
    Ignore IRDAD00014686.fits -> IRDAD00014686_cp.fits
    Ignore IRDAD00014687.fits -> IRDAD00014687_cp.fits
    Ignore IRDAD00014688.fits -> IRDAD00014688_cp.fits
    Ignore IRDAD00014689.fits -> IRDAD00014689_cp.fits
    Ignore IRDAD00014690.fits -> IRDAD00014690_cp.fits
    Ignore IRDAD00014691.fits -> IRDAD00014691_cp.fits
    Ignore IRDAD00014692.fits -> IRDAD00014692_cp.fits
    Ignore IRDAD00014693.fits -> IRDAD00014693_cp.fits
    Ignore IRDAD00014694.fits -> IRDAD00014694_cp.fits
    Ignore IRDAD00014695.fits -> IRDAD00014695_cp.fits
    Ignore IRDAD00014696.fits -> IRDAD00014696_cp.fits
    Ignore IRDAD00014697.fits -> IRDAD00014697_cp.fits
    Ignore IRDAD00014698.fits -> IRDAD00014698_cp.fits
    Ignore IRDAD00014699.fits -> IRDAD00014699_cp.fits
    Ignore IRDAD00014700.fits -> IRDAD00014700_cp.fits
    Ignore IRDAD00014701.fits -> IRDAD00014701_cp.fits
    Ignore IRDAD00014702.fits -> IRDAD00014702_cp.fits
    Ignore IRDAD00014703.fits -> IRDAD00014703_cp.fits
    Ignore IRDAD00014704.fits -> IRDAD00014704_cp.fits
    Ignore IRDAD00014705.fits -> IRDAD00014705_cp.fits
    Ignore IRDAD00014706.fits -> IRDAD00014706_cp.fits
    Ignore IRDAD00014707.fits -> IRDAD00014707_cp.fits
    Ignore IRDAD00014708.fits -> IRDAD00014708_cp.fits
    Ignore IRDAD00014709.fits -> IRDAD00014709_cp.fits
    Ignore IRDAD00014710.fits -> IRDAD00014710_cp.fits
    Ignore IRDAD00014711.fits -> IRDAD00014711_cp.fits
    Ignore IRDAD00014712.fits -> IRDAD00014712_cp.fits
    Ignore IRDAD00014713.fits -> IRDAD00014713_cp.fits
    Ignore IRDAD00014714.fits -> IRDAD00014714_cp.fits
    Ignore IRDAD00014715.fits -> IRDAD00014715_cp.fits
    Ignore IRDAD00014716.fits -> IRDAD00014716_cp.fits
    Ignore IRDAD00014717.fits -> IRDAD00014717_cp.fits
    Ignore IRDAD00014718.fits -> IRDAD00014718_cp.fits
    Ignore IRDAD00014719.fits -> IRDAD00014719_cp.fits
    Ignore IRDAD00014720.fits -> IRDAD00014720_cp.fits
    Ignore IRDAD00014721.fits -> IRDAD00014721_cp.fits
    Ignore IRDAD00014722.fits -> IRDAD00014722_cp.fits
    Ignore IRDAD00014723.fits -> IRDAD00014723_cp.fits
    Ignore IRDAD00014724.fits -> IRDAD00014724_cp.fits
    Ignore IRDAD00014725.fits -> IRDAD00014725_cp.fits
    Ignore IRDAD00014726.fits -> IRDAD00014726_cp.fits
    Ignore IRDAD00014727.fits -> IRDAD00014727_cp.fits
    Ignore IRDAD00014728.fits -> IRDAD00014728_cp.fits
    Ignore IRDAD00014729.fits -> IRDAD00014729_cp.fits
    Ignore IRDAD00014730.fits -> IRDAD00014730_cp.fits
    Ignore IRDAD00014731.fits -> IRDAD00014731_cp.fits
    Skipped 100 files because they already exists.


.. parsed-literal::

    0it [00:00, ?it/s]


.. parsed-literal::

    median combine:  _cp


.. parsed-literal::

    100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████| 100/100 [00:00<00:00, 890.20it/s]
    /Users/yuikasagi/miniforge3/envs/py39_pip/lib/python3.9/site-packages/numpy/lib/nanfunctions.py:1218: RuntimeWarning: All-NaN slice encountered
      r, k = function_base._ureduce(a, func=_nanmedian, axis=axis, out=out,


.. parsed-literal::

    fitsid: [41510]
    clean_pattern: output extension=_cp


.. parsed-literal::

    100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 21/21 [00:00<00:00, 102.09it/s]


.. parsed-literal::

    Ignore IRDA00041511.fits -> IRDA00041511_cp.fits


.. parsed-literal::

    0it [00:00, ?it/s]


Settings to display figures
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can set some parameters for plot(s). For example, you can get images
for several orders.

.. code:: ipython3

    import matplotlib
    matplotlib.use('tkagg')
    from pyird.utils.image_widget import image_1Dand2D
    import tkinter as tk
    
    ### SET PARAMETERS ###
    hotpix_mask = None # comment out this if you want to show hotpixels
    target.imcomb = False # set 'True' if you want to median combine images.
    wavcal_path = thar.anadir/('thar_%s_%s.fits'%(thar.band,thar.trace.mmf))
    
    ## additional parameters for plot
    vmin = -10
    vmax = 50
    scale = 'linear' # 'linear' or 'log'
    params = {'vmin':vmin,'vmax':vmax,'scale':scale}
    
    orders=[10,12] # orders to be plotted
    #######################
    
    ## Values needed for plotting
    rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot = target.flatten(extin='_cp',check=True,master_path=wavcal_path)


.. parsed-literal::

      0%|                                                                                                                           | 0/1 [00:00<?, ?it/s]
      0%|                                                                                                                          | 0/21 [00:00<?, ?it/s][A
     24%|███████████████████████████▏                                                                                      | 5/21 [00:00<00:00, 48.12it/s][A
     48%|█████████████████████████████████████████████████████▊                                                           | 10/21 [00:00<00:00, 44.28it/s][A
     71%|████████████████████████████████████████████████████████████████████████████████▋                                | 15/21 [00:00<00:00, 43.52it/s][A
    100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 21/21 [00:00<00:00, 43.50it/s]
      0%|                                                                                                                           | 0/1 [00:00<?, ?it/s]


Case 1. Plot Absorption Lines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``show_spec_to_image()`` will create two figures in one window for each
order.

-  The upper figure is the spectrum of a order, and the lower figure is
   the detector image.
-  You can zoom up both image.
-  When you press any key on the spectrum, corresponding positions on
   the detector will be plotted as ‘x’.

**Note**:

If you run in jupyter notebook, add ``%matplolib notebook`` and comment
out ``root.mainloop()``.

.. code:: ipython3

    ## show 1d spectrum and 2d image
    %matplotlib notebook
    for order in orders:
        print(order)
        ## draw window
        root = tk.Tk()
        window = image_1Dand2D(root,order=order,band=flat.band)
        window.show_spec_to_image(rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,wavcal_path=wavcal_path,hotpix_mask=hotpix_mask,**params)
    #root.mainloop()


.. parsed-literal::

    10



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <div id='a3c5d1a2-0f6f-4b78-bc61-258f84687329'></div>


.. parsed-literal::

    12



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <div id='9e2b2340-cfc1-4e73-83e5-61fd1687781f'></div>


Case 2. Plot Emission Lines
~~~~~~~~~~~~~~~~~~~~~~~~~~~

``show_emission_position()`` will be useful for the emission spectrum
(e.g. sky spectrum).

-  The upper figure is the detector image of one aperture, and the lower
   figure is the spectra of the order.
-  By fitting 2D gaussian to the emissions on the detector, the emission
   like signal and hotpixels are distinguished automatically.

**Note**:

If you run in jupyter notebook, add ``%matplolib notebook`` and comment
out ``root.mainloop()``.

.. code:: ipython3

    ## show positions of emissions on a detector image
    %matplotlib notebook
    for order in orders:
        ## draw window
        root2 = tk.Tk()
        window2 = image_1Dand2D(root2,order=order,band=flat.band)
        window2.show_emission_position(target,rsd,wav,mask,pixcoord,rotim,iys_plot,iye_plot,wavcal_path=wavcal_path,hotpix_mask=hotpix_mask,**params)
    #root2.mainloop()



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <div id='19180e1b-88e4-4265-b3a9-5acb65876303'></div>


.. parsed-literal::

    WARNING: The fit may be unsuccessful; check fit_info['message'] for more information. [astropy.modeling.fitting]



.. parsed-literal::

    <IPython.core.display.Javascript object>



.. raw:: html

    <div id='ba25dec4-f7a1-4eb1-ab08-533f4d1ccd8e'></div>


.. parsed-literal::

    WARNING: The fit may be unsuccessful; check fit_info['message'] for more information. [astropy.modeling.fitting]

