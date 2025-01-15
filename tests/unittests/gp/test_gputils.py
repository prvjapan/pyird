import pytest
import numpy as np
from pyird.gp.gputils import calc_coarsed_array 

@pytest.fixture
def sample_array():
    return np.array([
        [1., 2., 3., 4.],
        [5., 6., 7., 8.],
        [9., 10., 11., 12.],
        [13., 14., 15., 16.]
    ])

def test_calc_coarsed_array_basic(sample_array):
    """Test basic functionality without cube mode."""
    Ncor = 2
    result = calc_coarsed_array(sample_array, Ncor, cube=False)
    expected_result = np.array([
        [1.5, 3.5],
        [5.5, 7.5],
        [9.5, 11.5],
        [13.5, 15.5]
    ])
    np.testing.assert_almost_equal(result, expected_result, decimal=6)

def test_calc_coarsed_array_with_cube(sample_array):
    """Test functionality with cube mode enabled."""
    Ncor = 2
    result = calc_coarsed_array(sample_array, Ncor, cube=True)
    expected_result = np.array([
        [3.5, 5.5],
        [11.5, 13.5]
    ])
    np.testing.assert_almost_equal(result, expected_result, decimal=6)

def test_calc_coarsed_array_with_outliers():
    """Test array with outliers."""
    array = np.ones((100,100))
    array[1,1] = 1000.
    Ncor = 2
    result = calc_coarsed_array(array, Ncor, cube=True)
    expected_result = np.ones((50,50))
    np.testing.assert_almost_equal(result, expected_result, decimal=6)