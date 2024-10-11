import pytest
import numpy as np
from pyird.spec.rsdmat import multiorder_to_rsd, rsd_order_medfilt 

def test_multiorder_to_rsd():
    rawspec = [np.array([1.0, 2.0, 3.0]), np.array([4.0, 5.0, 6.0])]
    pixcoord = [np.array([0, 1, 2]), np.array([3, 4, 5])]
    npix = 10
    fill_value = np.nan
    expected_rsd = np.array([
        [1.0, np.nan],
        [2.0, np.nan],
        [3.0, np.nan],
        [np.nan, 4.0],
        [np.nan, 5.0],
        [np.nan, 6.0],
        [np.nan, np.nan],
        [np.nan, np.nan],
        [np.nan, np.nan],
        [np.nan, np.nan]
    ])
    rsd = multiorder_to_rsd(rawspec, pixcoord, npix=npix, fill_value=fill_value)
    np.testing.assert_array_equal(rsd, expected_rsd)

def test_rsd_order_medfilt():
    rsd = np.array([
        [1.0, 4.0],
        [2.0, 5.0],
        [3.0, 6.0],
        [7.0, 8.0],
        [9.0, 10.0]
    ])
    kernel_size = 3
    expected_filtered_rsd = np.array([
        [1.0, 4.0],
        [2.0, 5.0],
        [3.0, 6.0],
        [7.0, 8.0],
        [7.0, 8.0]
    ])  # Manually verified result of applying median filter
    rsd_filtered = rsd_order_medfilt(rsd, kernel_size=kernel_size)
    print(rsd_filtered)
    np.testing.assert_array_equal(rsd_filtered, expected_filtered_rsd)

if __name__ == '__main__':
    pytest.main()