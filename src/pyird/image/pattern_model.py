import numpy as np
from pyird.image.channel import image_to_channel_cube, channel_cube_to_image, eopixel_split, eopixel_combine
from pyird.plot.detector import show_profile

def median_XY_profile(calim, rm_nct=True, Ncor=64, show=True):
    """a simple readout pattern model.

    Note:
        This function assumes model = cmosub_med X + cmosub_med Y+ cmo, where cmos_med=cmo-subtracted median and cmo=channel median offset.

    Args:
        calim: masked image for read pattern calibration
        rm_nct: remove non-common trends of channel using GP
        Ncor: coarse graing number for rm_nct
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
        Nch=np.shape(model_channel_cube)[0]
        for i in range(0,Nch):
            cgda=model_channel_cube[i,:,:]
            calc_coarsed_array.coarse_gp(cgda, subarray, Ncor)
            
    model_image = channel_cube_to_image(model_channel_cube)

    return model_image
