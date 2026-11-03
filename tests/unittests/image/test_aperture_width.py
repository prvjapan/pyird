import numpy as np
import pytest

from pyird.image.aperture_width import resolve_aperture_width


@pytest.mark.parametrize(
    ("inst", "expected"),
    [
        ("IRD", (2, 4)),
        ("REACH", (2, 3)),
        ("IRCS", (4, 4)),
    ],
)
def test_resolve_aperture_width_uses_instrument_default(inst, expected):
    assert resolve_aperture_width(inst) == expected


def test_resolve_aperture_width_accepts_explicit_width():
    width = resolve_aperture_width(
        "UNKNOWN", [np.int64(3), np.int64(5)]
    )

    assert width == (3, 5)
    assert all(type(value) is int for value in width)


def test_resolve_aperture_width_rejects_unknown_instrument_without_width():
    with pytest.raises(ValueError, match="No default aperture width"):
        resolve_aperture_width("UNKNOWN")


@pytest.mark.parametrize(
    "width",
    [
        1,
        (1,),
        (1, 2, 3),
        (-1, 2),
        (1, -2),
        (1.5, 2),
        (True, 2),
        (1, False),
    ],
)
def test_resolve_aperture_width_rejects_invalid_width(width):
    with pytest.raises(ValueError, match="two non-negative integers"):
        resolve_aperture_width("IRD", width)
