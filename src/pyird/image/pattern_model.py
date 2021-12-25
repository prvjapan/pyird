
def median_XY_profile(calim):
    """
    
    Args:
        calim: masked image for read pattern calibration
    
    Returns:
        model pattern image

    """
    cal_channel_cube=image_to_channel_cube(calim,revert=True)
    cal_eotensor=eopixel_split(cal_channel_cube)
    _,nchan,xsize,ysize=np.shape(cal_eotensor)
    
    channel_median_offset=np.nanmedian(cal_eotensor,axis=(2,3))
    xprofile_offset_subtracted=np.nanmedian(cal_eotensor,axis=3)-channel_median_offset[:,:,np.newaxis]
    yprofile_offset_subtracted=np.nanmedian(cal_eotensor,axis=2)-channel_median_offset[:,:,np.newaxis]
    
    #median as model--------------------
    xprofile_offset_subtracted_model=np.nanmedian(xprofile_offset_subtracted,axis=1)
    yprofile_offset_subtracted_model=np.nanmedian(yprofile_offset_subtracted,axis=1)
    
    image_pattern_model=\
        xprofile_offset_subtracted_model[:,np.newaxis,:,np.newaxis]\
        +yprofile_offset_subtracted_model[:,np.newaxis,np.newaxis,:]\
        +channel_median_offset[:,:,np.newaxis,np.newaxis]
    return image_pattern_model
