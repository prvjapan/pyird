Handling Hot Pixels and Order Edges
====================================

This page explains how hot pixels and order edges are treated in PyIRD.

Hot Pixels → NaN or Interpolation
----------------------------------

The hotpix_mask can be created in two ways:

1. Generated from DARK frames
2. Loaded from a prepared mask (see :doc:`./other_settings`)

This mask is applied in the following functions of the Stream2D class:

- ``clean_pattern()``: Applies the hot pixel mask when modeling detector patterns.
- ``flatten()``: Replaces hot pixels with NaN or interpolated values using a spline function.
- ``apnormalize()``: Same as flatten().
- ``apext_flatfield()``: Same as flatten().

Note that the hot pixel mask does not affect the extracted spectrum in ``clean_pattern()``, while in the other three functions it can affect the results through interpolation via ``pyird.image.hotpix.apply_hotpixel_orderedge_mask``.

For the last three functions, you can set the variable ``hotpix_mode``.
This variable controls how hot pixels are treated:

- If ``hotpix_mode="nan"``, bad pixels are replaced with ``np.nan``.
- If ``hotpix_mode="interp"``, they are replaced using spline interpolation.

.. note::

    ``hotpix_mode`` is available in PyIRD version > 1.2, and the default is ``"nan"``.
    In previous versions, the code used ``hotpix_mode="interp"``.

Check the Hot Pixel Masked Positions
"""""""""""""""""""""""""""""""""""""
If you want to check which pixels are masked by hotpix_mask, you can extract their spectra by inserting the following code into IRD_stream.py.

.. code:: python

    # Example
    # Insert at the end of IRD_stream.py
    dark.trace = trace_mmf
    output_fl = f"hotpix_mask_fl_{band}_m{mmf[-1]}.fits"
    dark.flatten_im(image=hotpix_mask, trace_path=None, output_path=output_fl)
    output_w = f"whotpix_mask_{band}_m{mmf[-1]}.dat"
    dark.dispcor_im(input_path=output_fl, master_path=thar.anadir, output_path=output_w)

Order Edges → Zero
-------------------
Flux at order edges, where the aperture trace fails, is automatically set to zero.
