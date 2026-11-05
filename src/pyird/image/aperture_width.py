"""Utilities for resolving aperture extraction widths."""

from numbers import Integral


DEFAULT_APERTURE_WIDTH = {
    "IRD": (2, 4),
    "REACH": (2, 3),
    "IRCS": (4, 4),
}


def resolve_aperture_width(inst, width=None):
    """Return a validated ``(width_start, width_end)`` tuple."""
    if width is None:
        try:
            return DEFAULT_APERTURE_WIDTH[inst]
        except KeyError as exc:
            raise ValueError(
                f"No default aperture width is defined for inst={inst!r}."
            ) from exc

    try:
        width_start, width_end = width
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "width must contain two non-negative integers."
        ) from exc

    resolved = (width_start, width_end)
    if any(
        isinstance(value, bool) or not isinstance(value, Integral) or value < 0
        for value in resolved
    ):
        raise ValueError("width must contain two non-negative integers.")

    return tuple(int(value) for value in resolved)
