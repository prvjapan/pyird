import numpy as np
from astropy.stats import sigma_clip
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.plot.detector import show_profile

def cal_nct(nctrend_im,margin_npixel=4,Ncor=64, sigma=0.1, xscale=32, yscale=64, cube=False):
    """
    Calculate the non-common trend model using Gaussian Process.

    Args:
        nctrend_im (ndarray): Input image for calculating non-common trend.
        margin_npixel (int): Margin pixels to be excluded.
        Ncor (int): Coarsening factor.
        sigma (float): Sigma parameter for GP.
        xscale (int): Scale factor for X-axis.
        yscale (int): Scale factor for Y-axis.
        cube (bool): If True, coarsen with cube structure.

    Returns:
        nctrend_model (ndarray): Non-common trend model.
    """
    from pyird.gp.gputils import calc_coarsed_array
    from gpkron.gp2d import GP2D, RBF
    subarray = np.zeros_like(nctrend_im)
    subarray = subarray[:, margin_npixel:-margin_npixel]

    coarsed_array = calc_coarsed_array(nctrend_im, Ncor, cube=cube)
    coarsed_array[coarsed_array !=coarsed_array] = np.nanmedian(coarsed_array)

    nctrend_model = GP2D(coarsed_array, RBF, sigma, (xscale, yscale), pshape=np.shape(subarray))
    return nctrend_model

def median_profile_xy(cal_eotensor_filtered, channel_median_offset, channel_axis=1):
    """
    Calculate the median profile of all even/odd channels while preserving the characteristics of each pixel.
    
    Parameters:
        cal_eotensor_filtered: numpy array, input tensor.
        channel_median_offset: numpy array, offset to subtract from each channel.
        channel_axis: int, axis representing the channels (default is 1).
    
    Returns:
        offset_subtracted: numpy array, input tensor with channel median offset subtracted.
        xyprofile_offset_subtracted_model: numpy array, median profile model preserving pixel profiles.
    """
    # Subtract channel median offset
    offset_subtracted = cal_eotensor_filtered - channel_median_offset[:, :, np.newaxis, np.newaxis]

    # Calculate median profiles for even and odd channels
    median_even_channels = np.nanmedian(offset_subtracted[:, ::2, :, :], axis=channel_axis)
    median_odd_channels = np.nanmedian(offset_subtracted[:, 1::2, :, :], axis=channel_axis)

    # Initialize the model tensor with the same shape as the input tensor
    xyprofile_offset_subtracted_model = np.zeros_like(cal_eotensor_filtered)

    # Fill in the model tensor for even and odd channels
    for i in range(2):
        xyprofile_offset_subtracted_model[i::2, ::2, :, :] = median_even_channels[i, np.newaxis, :, :]
        xyprofile_offset_subtracted_model[i::2, 1::2, :, :] = median_odd_channels[i, np.newaxis, :, :]

    return offset_subtracted, xyprofile_offset_subtracted_model

def pad_nans_with_xy_median(offset_subtracted, channel_median_offset, image_pattern_model_eotensor, channel_axis=1, x_axis=2, y_axis=3):
    """
    Pad NaN values in the input tensor with the median profile calculated from x and y axes.
    
    Parameters:
        offset_subtracted: numpy array, input tensor with subtracted offset.
        channel_median_offset: numpy array, offset to add back to the channels.
        image_pattern_model_eotensor: numpy array, tensor with NaNs to be filled.
        channel_axis: int, axis representing the channels (default is 1).
        x_axis: int, axis representing the x-axis (default is 2).
        y_axis: int, axis representing the y-axis (default is 3).
    
    Returns:
        image_pattern_model_eotensor: numpy array, tensor with NaNs filled.
    """
    # Calculate median profiles along x and y axes
    xprofile_offset_subtracted = np.nanmedian(offset_subtracted, axis=y_axis)
    yprofile_offset_subtracted = np.nanmedian(offset_subtracted, axis=x_axis)

    # Calculate median profiles for each channel
    xprofile_offset_subtracted_model = np.nanmedian(xprofile_offset_subtracted, axis=channel_axis)
    yprofile_offset_subtracted_model = np.nanmedian(yprofile_offset_subtracted, axis=channel_axis)

    # Construct the grid to fill NaNs
    xygrid_offset_subtracted_model = (
        xprofile_offset_subtracted_model[:, np.newaxis, :, np.newaxis]
        + yprofile_offset_subtracted_model[:, np.newaxis, np.newaxis, :]
        + channel_median_offset[:, :, np.newaxis, np.newaxis]
    )

    # Fill NaN values in the image pattern model tensor
    nan_mask = np.isnan(image_pattern_model_eotensor)
    image_pattern_model_eotensor[nan_mask] = xygrid_offset_subtracted_model[nan_mask]
    return image_pattern_model_eotensor

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
        nctrend_model_sc = cal_nct(calim0, Ncor=64, xscale=512, yscale=512, cube=True, margin_npixel=margin_npixel)
        sc_model[:, margin_npixel:-margin_npixel] = nctrend_model_sc
        model_image = sc_model
        calim = calim0 - model_image
    else:
        calim = calim0

    ## channel cube and e/o tensor
    cal_channel_cube = image_to_channel_cube(calim, revert=True)
    cal_eotensor = eopixel_split(cal_channel_cube) ## eo tensor (Neo x Nchannel x Nx x Ny = 2 x 32 x 32 x 2048)
    eo_axis, channel_axis, column, row = 0, 1, 2, 3

    ## reject outliers
    cal_eotensor_cut = cal_eotensor[:, :, :, margin_npixel:-margin_npixel]
    mask = sigma_clip(cal_eotensor_cut,maxiters=1).mask
    cal_eotensor_filtered = cal_eotensor.copy()
    cal_eotensor_filtered[:, :, :, margin_npixel:-margin_npixel][mask] = np.nan

    ## overall median profile (offset for each channel)
    channel_median_offset = np.nanmedian(cal_eotensor_filtered, axis=(column, row)) # -> Neo x Nchannel

    ## median of all even/odd channels for the offset subtracted image
    offset_subtracted, xyprofile_offset_subtracted_model = median_profile_xy(cal_eotensor_filtered, channel_median_offset, channel_axis=channel_axis)
    ## median profile = channel offset + median pixel profile
    image_pattern_model_eotensor = xyprofile_offset_subtracted_model + channel_median_offset[:, :, np.newaxis, np.newaxis]

    ## pad nans (at aperture mask) with {median profile of x pixels + median profile of y pixels} 
    if np.sum(np.isnan(image_pattern_model_eotensor)):
        image_pattern_model_eotensor = pad_nans_with_xy_median(offset_subtracted, channel_median_offset, image_pattern_model_eotensor, 
                                                               channel_axis=channel_axis, x_axis=column, y_axis=row)
    
    model_channel_cube = eopixel_combine(image_pattern_model_eotensor)

    if rm_nct:
        Nch = np.shape(model_channel_cube)[0]
        for i in range(0, Nch):
            nctrend = cal_channel_cube[i, :, :] - model_channel_cube[i, :, :]

            nctrend_model = cal_nct(nctrend, margin_npixel=margin_npixel)
            model_channel_cube[i, :, margin_npixel:-margin_npixel] = model_channel_cube[i, :, margin_npixel:-margin_npixel] + nctrend_model

        model_image += channel_cube_to_image(model_channel_cube)
    else:
        model_image = channel_cube_to_image(model_channel_cube)

    ## stripe model
    npix = calim0.shape[0]
    corrected_im_tmp = calim0 - model_image
    stripe = np.nanmedian(corrected_im_tmp, axis=0) ## column direction
    stripe = np.array([stripe]*npix)
    model_image += stripe

    return model_image
