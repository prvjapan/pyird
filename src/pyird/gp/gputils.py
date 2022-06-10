import numpy as np

def calc_coarsed_array(array, Ncor):
    """coarsed array."""
    each = (array-np.median(array))
    mask = np.abs(each) > 5.0*np.std(each)
    marray = np.copy(array)
    marray[mask] = None
    stacked_array = []
    for j in range(0, Ncor):
        stacked_array.append(marray[:, j::Ncor])
    stacked_array = np.array(stacked_array)
    coarsed_array = np.nanmedian(stacked_array, axis=0)
    return coarsed_array


if __name__ == '__main__':

    #calc_coarsed_array(array, subarray, xscale, yscale, Ncor)
    #coarsed_array[coarsed_array != coarsed_array] = np.nanmedian(coarsed_array)
    sigma = 0.001
    #GLP = gp2d.GP2Dcross(coarsed_array, subarray, sigma, xscale, yscale)
