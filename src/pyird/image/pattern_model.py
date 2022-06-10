import numpy as np
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.plot.detector import show_profile


def median_XY_profile(calim, rm_nct=True, Ncor=64, margin_npixel=4, sigma=0.1, xscale=32, yscale=64, show=True):
    """a simple readout pattern model.

    Note:
        This function assumes model = cmosub_med X + cmosub_med Y+ cmo, where cmos_med=cmo-subtracted median and cmo=channel median offset. When rm_cnt=True option corrects non-common trends of channels using a 2D GP. See #10 (https://github.com/prvjapan/pyird/issues/10).

    Args:
        calim: masked image for read pattern calibration
        rm_nct: remove non-common trends of channel using a GP
        Ncor: coarse graing number for rm_nct
        margin_npixel: # of pixels for detector margin
        sigma: GP diagonal component
        xscale: GP x scale
        yscale: GP y scale
        show: showing profile

    Returns:
        model pattern image
    """
    cal_channel_cube = image_to_channel_cube(calim, revert=True)
    cal_eotensor = eopixel_split(cal_channel_cube)
    _, nchan, xsize, ysize = np.shape(cal_eotensor)

    channel_median_offset = np.nanmedian(cal_eotensor, axis=(2, 3))
    xprofile_offset_subtracted = np.nanmedian(
        cal_eotensor, axis=3)-channel_median_offset[:, :, np.newaxis]
    yprofile_offset_subtracted = np.nanmedian(
        cal_eotensor, axis=2)-channel_median_offset[:, :, np.newaxis]

    # median as model--------------------
    xprofile_offset_subtracted_model = np.nanmedian(
        xprofile_offset_subtracted, axis=1)
    yprofile_offset_subtracted_model = np.nanmedian(
        yprofile_offset_subtracted, axis=1)

    if show:
        show_profile(xprofile_offset_subtracted, yprofile_offset_subtracted,
                     xprofile_offset_subtracted_model, yprofile_offset_subtracted_model)

    image_pattern_model_eotensor =\
        xprofile_offset_subtracted_model[:, np.newaxis, :, np.newaxis]\
        + yprofile_offset_subtracted_model[:, np.newaxis, np.newaxis, :]\
        + channel_median_offset[:, :, np.newaxis, np.newaxis]

    model_channel_cube = eopixel_combine(image_pattern_model_eotensor)
    if rm_nct:
        from pyird.gp.gputils import calc_coarsed_array
        from gpkron.gp2d import GP2D, RBF

        Nch = np.shape(model_channel_cube)[0]
        for i in range(0, Nch):
            nctrend = cal_channel_cube[i, :, :]-model_channel_cube[i, :, :]

            subarray = np.zeros_like(nctrend)
            subarray = subarray[:, margin_npixel:-margin_npixel]

            coarsed_array = calc_coarsed_array(nctrend, Ncor)
            coarsed_array[coarsed_array !=
                          coarsed_array] = np.nanmedian(coarsed_array)

            nctrend_model = GP2D(coarsed_array, RBF, sigma,
                                 (xscale, yscale), pshape=np.shape(subarray))
            model_channel_cube[i, :, margin_npixel:-
                               margin_npixel] = model_channel_cube[i, :, margin_npixel:-margin_npixel]+nctrend_model

    model_image = channel_cube_to_image(model_channel_cube)

    return model_image
