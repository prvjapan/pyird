"""Read and write IRAF APEXTRACT aperture database files."""
from collections.abc import Iterable
from datetime import datetime
from os import PathLike
from pathlib import Path

import numpy as np

__all__ = ["read_trace_file", "read_trace_file_one", "write_trace_file"]


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
    if isinstance(filelist, (str, bytes, PathLike)):
        filenames = [filelist]
    elif isinstance(filelist, Iterable):
        filenames = list(filelist)
    else:
        raise TypeError("filelist must be a path or an iterable of paths")

    if not filenames:
        raise ValueError("filelist must contain at least one trace file")

    y0, interp_function, xmin, xmax, coeff = read_trace_file_one(
        filenames[0], finalize=False)
    for filename in filenames[1:]:
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

    with open(filename, encoding="utf-8") as f:
        cont = f.readlines()
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


def write_trace_file(trace, stream, background="none", rot180=False):
    """Write a trace to an IRAF APEXTRACT database file.

    Parameters
    ----------
    trace : TraceAperture
        Aperture trace to serialize.
    stream : Stream2D
        Supplies ``anadir`` and ``streamid`` for the output path.
    background : {"none", "default"}, optional
        Whether to include the standard PyIRD background definition.
    rot180 : bool, optional
        Transform coordinates back from a 180-degree rotated image.
    """
    if background not in ("none", "default"):
        raise ValueError("background must be either 'none' or 'default'.")

    n_disp, n_spatial = map(int, trace.mask_shape)
    if n_disp <= 0 or n_spatial <= 0:
        raise ValueError("mask_shape must contain positive integers.")
    if getattr(trace.trace_function, "__name__", None) != "trace_legendre":
        raise NotImplementedError(
            "write_trace_file currently supports only trace_legendre."
        )

    y0 = np.asarray(trace.y0, dtype=float)
    xmin = np.asarray(trace.xmin, dtype=float)
    xmax = np.asarray(trace.xmax, dtype=float)
    coeff = [np.asarray(c, dtype=float) for c in trace.coeff]
    nap = len(y0)
    if not (len(xmin) == nap and len(xmax) == nap and len(coeff) == nap):
        raise ValueError(
            "y0, xmin, xmax, and coeff must have the same number of apertures."
        )

    if trace.width is None:
        default_widths = {"IRD": (2, 4), "REACH": (2, 3)}
        try:
            width_str, width_end = default_widths[trace.inst]
        except KeyError as exc:
            raise ValueError(
                "No default aperture width is defined for "
                f"inst='{trace.inst}'. Please set trace.width."
            ) from exc
    else:
        if len(trace.width) != 2:
            raise ValueError("TraceAperture.width must be [width_str, width_end].")
        width_str, width_end = trace.width

    if rot180:
        aperture_low, aperture_high = -float(width_end), float(width_str)
        y0_iraf = n_spatial + 1.0 - y0
        xmin_iraf, xmax_iraf = n_disp - xmax, n_disp - xmin
        coeff_iraf = []
        for values in coeff:
            degree = np.arange(len(values))
            coeff_iraf.append(values * (-1.0) ** (degree + 1))
    else:
        aperture_low, aperture_high = -float(width_str), float(width_end)
        y0_iraf = y0.copy()
        xmin_iraf, xmax_iraf = xmin + 1.0, xmax + 1.0
        coeff_iraf = [values.copy() for values in coeff]

    aperture_order = np.argsort(y0_iraf, kind="stable")
    disp_center = (n_disp + 1) // 2
    disp_low, disp_high = 1 - disp_center, n_disp - disp_center
    image = f"{stream.streamid}_{trace.mmf}"
    database_dir = Path(stream.anadir) / "database"
    database_dir.mkdir(parents=True, exist_ok=True)
    output_path = database_dir / f"ap{image}"

    def _fmt(value):
        return f"{float(value):.12g}"

    with open(output_path, "w", encoding="utf-8") as handle:
        for aperture, index in enumerate(aperture_order, start=1):
            center = y0_iraf[index]
            values = coeff_iraf[index]
            timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
            handle.write(f"# {timestamp}\n")
            handle.write(
                f"begin\taperture {image} {aperture} "
                f"{_fmt(center)} {_fmt(disp_center)}\n"
            )
            handle.write(f"\timage\t{image}\n")
            handle.write(f"\taperture\t{aperture}\n")
            handle.write(f"\tbeam\t{aperture}\n")
            handle.write(f"\tcenter\t{_fmt(center)} {_fmt(disp_center)}\n")
            handle.write(f"\tlow\t{_fmt(aperture_low)} {_fmt(disp_low)}\n")
            handle.write(f"\thigh\t{_fmt(aperture_high)} {_fmt(disp_high)}\n")
            if background == "default":
                handle.write("\tbackground\n")
                handle.write("\t\txmin -10.\n\t\txmax 10.\n")
                handle.write("\t\tfunction chebyshev\n\t\torder 1\n")
                handle.write("\t\tsample -10:-6,6:10\n\t\tnaverage -3\n")
                handle.write("\t\tniterate 0\n\t\tlow_reject 3.\n")
                handle.write("\t\thigh_reject 3.\n\t\tgrow 0.\n")
            handle.write("\taxis\t1\n")
            handle.write(f"\tcurve\t{4 + len(values)}\n")
            handle.write("\t\t2.\n")
            handle.write(f"\t\t{len(values)}.\n")
            handle.write(f"\t\t{_fmt(xmin_iraf[index])}\n")
            handle.write(f"\t\t{_fmt(xmax_iraf[index])}\n")
            for value in values:
                handle.write(f"\t\t{_fmt(value)}\n")
            handle.write("\n")

    if trace.info:
        print(f"Created IRAF aperture database: {output_path}")
    return output_path


if __name__ == '__main__':
    import importlib
    pathA = (importlib.resources.files('pyird').joinpath('data/samples/aprefA'))
    pathC = (importlib.resources.files('pyird').joinpath('data/samples/aprefC'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file([pathA, pathC])
    print(len(y0))
