from types import SimpleNamespace

import numpy as np

from pyird.image import pattern_model
from pyird.utils import irdstream


def test_clean_pattern_rotates_y_mask_for_every_frame(monkeypatch, tmp_path):
    mask = np.zeros((2, 3), dtype=bool)
    mask[0, 0] = True
    seen_masks = []
    stream = SimpleNamespace(
        band="y",
        extension="",
        rawdir=tmp_path,
        anadir=tmp_path,
        print_if_info_is_true=lambda _: None,
        check_existence=lambda *_: (["a.fits", "b.fits"], ["a_cp.fits", "b_cp.fits"]),
        detector_handling=lambda image, mode: image,
    )

    def record_mask(image):
        seen_masks.append(np.isnan(image))
        return np.zeros_like(image)

    monkeypatch.setattr(
        irdstream, "load_fits_data_header", lambda _: (np.ones(mask.shape), {})
    )
    monkeypatch.setattr(irdstream, "write_fits_data_header", lambda *_: None)
    monkeypatch.setattr(pattern_model, "median_XY_profile", record_mask)

    irdstream.Stream2D.clean_pattern(stream, trace_mask=mask)

    expected = mask[::-1, ::-1]
    np.testing.assert_array_equal(seen_masks, [expected, expected])
