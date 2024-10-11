import pytest
import numpy as np
from pyird.image.trace_function import trace_legendre


def test_trace_legendre():
    # Test inputs
    x = [np.array([1, 2, 3]), np.array([4, 5, 6])]
    y0 = [1, 2]
    xmin = [0, 3]
    xmax = [5, 7]
    coeff = [[0.5, -0.2, 0.3], [1.0, 0.0, -0.5]]

    # Expected output (calculated manually or based on the expected behavior of the function)
    expected_trace_lines = [
        np.array([0.632, 0.408, 0.328]),
        np.array([2.0625, 2.25  , 2.0625])
    ]

    # Run the function
    result = trace_legendre(x, y0, xmin, xmax, coeff)

    # Assert each result against expected output
    for res, exp in zip(result, expected_trace_lines):
        np.testing.assert_almost_equal(res, exp, decimal=2)


if __name__ == "__main__":
    pytest.main()