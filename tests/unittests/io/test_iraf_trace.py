import pytest
import importlib
from types import SimpleNamespace

from pyird.image.trace_function import trace_legendre
from pyird.io.iraf_trace import read_trace_file, write_trace_file
from pyird.utils.aperture import TraceAperture
import numpy as np

def test_read_trace_file():
    path = (importlib.resources.files('pyird').joinpath('data/samples/aprefA'))
    y0, interp_function, xmin, xmax, coeff = read_trace_file(path)
    assert np.sum(y0) == 13517.2013


def test_read_trace_file_accepts_an_iterable():
    path = importlib.resources.files("pyird").joinpath("data/samples/aprefA")
    expected = read_trace_file(path)
    actual = read_trace_file(iter([path]))

    for actual_value, expected_value in zip(actual[:4], expected[:4]):
        np.testing.assert_equal(actual_value, expected_value)
    assert actual[4] == expected[4]


def test_read_trace_file_rejects_empty_iterable():
    with pytest.raises(ValueError, match="at least one"):
        read_trace_file([])


def test_write_trace_file_round_trip(tmp_path):
    aperture = TraceAperture(
        trace_legendre,
        y0=[10.0, 20.0],
        xmin=[0, 1],
        xmax=[99, 98],
        coeff=[[1.0, 0.25], [2.0, -0.5]],
        inst="IRD",
        mask_shape=(100, 50),
    )
    stream = SimpleNamespace(anadir=tmp_path, streamid="test")

    output_path = write_trace_file(aperture, stream, background="none")
    y0, function, xmin, xmax, coeff = read_trace_file(output_path)

    assert output_path == tmp_path / "database" / "aptest_mmf12"
    np.testing.assert_equal(y0, [10.0, 20.0])
    np.testing.assert_equal(function, [2, 2])
    np.testing.assert_equal(xmin, [1, 2])
    np.testing.assert_equal(xmax, [100, 99])
    np.testing.assert_equal(coeff, [[1.0, 0.25], [2.0, -0.5]])


def test_trace_aperture_write_iraf_delegates_to_writer(tmp_path):
    aperture = TraceAperture(
        trace_legendre,
        y0=[10.0],
        xmin=[0],
        xmax=[99],
        coeff=[[1.0]],
        inst="IRD",
        mask_shape=(100, 50),
    )
    stream = SimpleNamespace(anadir=tmp_path, streamid="test")

    output_path = aperture.write_iraf(stream, background="none")

    assert output_path.is_file()

if __name__ == '__main__':
    test_read_trace_file()
