import numpy as np
from astropy.stats import sigma_clip
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.plot.detector import show_profile

def cal_nct(nctrend_im,margin_npixel=4,Ncor=64, sigma=0.1, xscale=32, yscale=64, cube=False):
    from pyird.gp.gputils import calc_coarsed_array
    from gpkron.gp2d import GP2D, RBF
    subarray = np.zeros_like(nctrend_im)
    subarray = subarray[:, margin_npixel:-margin_npixel]

    coarsed_array = calc_coarsed_array(nctrend_im, Ncor, cube=cube)
    coarsed_array[coarsed_array !=coarsed_array] = np.nanmedian(coarsed_array)

    nctrend_model = GP2D(coarsed_array, RBF, sigma, (xscale, yscale), pshape=np.shape(subarray))
    return nctrend_model

def median_XY_profile(calim0, rm_nct=True, margin_npixel=4):
    """a simple readout pattern model. (revised by Y.K., May 2023)

    Note:
        This function assumes model = cmosub_med X-Y + cmo, where cmos_med=cmo-subtracted median and cmo=channel median offset. When rm_cnt=True option corrects non-common trends of channels using a 2D GP. See #10 (https://github.com/prvjapan/pyird/issues/10).

    Args:
        calim0: masked image for read pattern calibration
        rm_nct: remove non-common trends of channel using a GP
        margin_npixel: # of pixels for detector margin

   Returns:
        model pattern image
    """
    import warnings
    warnings.simplefilter('ignore')

    ## overall trend model (scatter light)
    if rm_nct:
        sc_model = calim0.copy()
        nctrend_model_sc = cal_nct(calim0,Ncor=64,xscale=512,yscale=512,cube=True)
        sc_model[:,margin_npixel:-margin_npixel] = nctrend_model_sc
        model_image = sc_model
        calim = calim0 - model_image
    else:
        calim = calim0

    ## channel cube and e/o tensor
    cal_channel_cube = image_to_channel_cube(calim, revert=True)
    cal_eotensor = eopixel_split(cal_channel_cube) ## eo tensor (2 [e/o] x Nchannel x xsize x ysize)
    nchan, xsize, ysize = np.shape(cal_channel_cube)

    ## reject outliers
    cal_eotensor_cut = cal_eotensor[:,:,:,margin_npixel:-margin_npixel]
    mask = sigma_clip(cal_eotensor_cut,maxiters=1).mask
    cal_eotensor_filtered = cal_eotensor.copy()
    cal_eotensor_filtered[:,:,:,margin_npixel:-margin_npixel][mask] = np.nan

    ## cmo (channel median offset)
    channel_median_offset = np.nanmedian(cal_eotensor_filtered, axis=(2, 3))

    ## cmosub_med X-Y
    offset_subtracted = cal_eotensor_filtered - channel_median_offset[:, :, np.newaxis, np.newaxis]
    xyprofile_offset_subtracted_model_echan = np.nanmedian(offset_subtracted[:,::2,:,:], axis=1)
    xyprofile_offset_subtracted_model_ochan = np.nanmedian(offset_subtracted[:,1::2,:,:], axis=1)
    xyprofile_offset_subtracted_model = np.zeros(cal_eotensor_filtered.shape)
    xyprofile_offset_subtracted_model[::2,::2,:,:] = np.array([xyprofile_offset_subtracted_model_echan[0,:,:]]*16).reshape(1,16,32,2048)
    xyprofile_offset_subtracted_model[1::2,::2,:,:] = np.array([xyprofile_offset_subtracted_model_echan[1,:,:]]*16).reshape(1,16,32,2048)
    xyprofile_offset_subtracted_model[::2,1::2,:,:] = np.array([xyprofile_offset_subtracted_model_ochan[0,:,:]]*16).reshape(1,16,32,2048)
    xyprofile_offset_subtracted_model[1::2,1::2,:,:] = np.array([xyprofile_offset_subtracted_model_ochan[1,:,:]]*16).reshape(1,16,32,2048)

    image_pattern_model_eotensor = xyprofile_offset_subtracted_model + channel_median_offset[:, :, np.newaxis, np.newaxis]

    ## pad nans with {cmosub_med X + cmosub_med Y}
    xprofile_offset_subtracted = np.nanmedian(offset_subtracted, axis=3)
    yprofile_offset_subtracted = np.nanmedian(offset_subtracted, axis=2)
    xprofile_offset_subtracted_model = np.nanmedian(xprofile_offset_subtracted, axis=1)
    yprofile_offset_subtracted_model = np.nanmedian(yprofile_offset_subtracted, axis=1)
    xygrid_offset_subtracted_model =\
        xprofile_offset_subtracted_model[:, np.newaxis, :, np.newaxis]\
        + yprofile_offset_subtracted_model[:, np.newaxis, np.newaxis, :]\
        + channel_median_offset[:, :, np.newaxis, np.newaxis]
    image_pattern_model_eotensor[np.isnan(image_pattern_model_eotensor)] =\
        xygrid_offset_subtracted_model[np.isnan(image_pattern_model_eotensor)]

    model_channel_cube = eopixel_combine(image_pattern_model_eotensor)

    if rm_nct:
        Nch = np.shape(model_channel_cube)[0]
        for i in range(0, Nch):
            nctrend = cal_channel_cube[i, :, :]-model_channel_cube[i, :, :]

            nctrend_model = cal_nct(nctrend)
            model_channel_cube[i, :, margin_npixel:-
                               margin_npixel] = model_channel_cube[i, :, margin_npixel:-margin_npixel]+nctrend_model

        model_image += channel_cube_to_image(model_channel_cube)
    else:
        model_image = channel_cube_to_image(model_channel_cube)

    ## stripe model
    corrected_im_tmp = calim0 - model_image
    stripe = np.nanmedian(corrected_im_tmp,axis=0) ## 0
    stripe = np.array([stripe]*2048)
    model_image += stripe

    return model_image
