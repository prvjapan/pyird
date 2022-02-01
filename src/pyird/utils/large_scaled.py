from pyird import gp2d
import numpy as np
####


def check_mask_filling(mask, Ncor=64):
    stacked_array = []
    for i in range(0, Ncor):
        for j in range(0, Ncor):
            stacked_array.append(mask[i::Ncor, j::Ncor])
    stacked_array = np.array(stacked_array)
    course_mask = np.sum(stacked_array, axis=0)
    if len(course_mask[course_mask == Ncor*Ncor]) > 0:
        print('Warning: Too many mask. Fill them by smoothing.')
        check = False
    else:
        check = True
    return check


def get_LSD(img, gpscale=1024, Ncor=64, sigma=0.001):
    stacked_array = []
    for i in range(0, Ncor):
        for j in range(0, Ncor):
            stacked_array.append(img[i::Ncor, j::Ncor])
    stacked_array = np.array(stacked_array)
    course_img = np.nanmedian(stacked_array, axis=0)
    LSD = gp2d.GP2Dcross(course_img, img, sigma, gpscale, gpscale)
    return LSD
