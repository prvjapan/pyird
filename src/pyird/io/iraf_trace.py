"""Load IRAF-made aperture file."""
import numpy as np


def finalize_trace(interp_function, xmin, xmax):
    interp_function = np.array(interp_function, dtype=int)
    xmin = np.array(xmin, dtype=int)
    xmax = np.array(xmax, dtype=int)
    return interp_function, xmin, xmax


def read_trace_file(filelist):
    """
    Args:
       filelist: trace files list

    Returns:
        y0, interp_function, xmin, xmax, coeff

    Examples:        
        >>> pathC=(pkg_resources.resource_filename('pyird', "data/samples/aprefC"))
        >>> path_c=(pkg_resources.resource_filename('pyird', "data/samples/apref_c"))
        >>> y0, interp_function, xmin, xmax, coeff=read_trace_file([pathC,path_c])


    """
    try:
        return read_trace_file_one(filelist)
    except:
        nf = len(filelist)

    if nf == 1:
        return read_trace_file_one(filelist[0])

    y0, interp_function, xmin, xmax, coeff = read_trace_file_one(
        filelist[0], finalize=False)
    for filename in filelist[1:]:
        y0, interp_function, xmin, xmax, coeff = read_trace_file_one(filename, finalize=False,
                                                                     y0=y0,
                                                                     interp_function=interp_function,
                                                                     xmin=xmin,
                                                                     xmax=xmax,
                                                                     coeff=coeff)
    interp_function, xmin, xmax = finalize_trace(interp_function, xmin, xmax)

    return y0, interp_function, xmin, xmax, coeff


def read_trace_file_one(filename, finalize=True, y0=None, interp_function=None, xmin=None, xmax=None, coeff=None):
    """
    Args:
       filename:trace file
       finalize: if you do not read trace file anymore, specify True, otherwise False.

    Returns:
        y0, interp_function, xmin, xmax, coeff

    """

    with open(filename) as f:
        cont = f.readlines()
    f.close()
    norder = 0
    read_curve = False
    curvepar = []

    if y0 is None:
        y0 = []
        interp_function = []
        xmin = []
        xmax = []
        coeff = []

    for line in cont:
        arr = line.split()
        if len(arr) > 0:
            # identify begin
            if arr[0] == 'begin':
                # reset read_curve
                if len(curvepar) > 0:
                    interp_function.append(float(curvepar[0]))
                    nlegendreorder = int(curvepar[1])
                    xmin.append(float(curvepar[2]))
                    xmax.append(float(curvepar[3]))
                    coeffval = [float(s) for s in curvepar[4:4+nlegendreorder]]
                    coeff.append(coeffval)
                    read_curve = False

                norder = norder+1
                y0.append(float(arr[4]))
            elif arr[0] == 'curve':
                read_curve = True
                curvepar = []
            elif read_curve and arr[0] != '#':
                curvepar.append(float(arr[0]))

    interp_function.append(float(curvepar[0]))
    xmin.append(float(curvepar[2]))
    xmax.append(float(curvepar[3]))
    nlegendreorder = int(curvepar[1])
    coeffval = [float(s) for s in curvepar[4:4+nlegendreorder]]
    coeff.append(coeffval)

    if finalize == True:
        interp_function, xmin, xmax = finalize_trace(
            interp_function, xmin, xmax)

    return y0, interp_function, xmin, xmax, coeff


if __name__ == '__main__':
    import pkg_resources
    pathA = (pkg_resources.resource_filename('pyird', 'data/samples/aprefA'))
    pathC = (pkg_resources.resource_filename('pyird', 'data/samples/aprefC'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file([pathA, pathC])
    print(len(y0))
