import numpy as np
import pytest
from pyird.image.pattern_model import median_XY_profile 

@pytest.fixture
def create_test_data():
    np.random.seed(42)
    # Create synthetic data for testing
    calim0 = np.random.normal(loc=100., scale=10., size=(2048, 2048))
    # Introduce NaNs to simulate a realistic input
    calim0[100:-100, 1500:1505] = np.nan
    import matplotlib.pyplot as plt
    plt.imshow(calim0, origin="lower")
    plt.show()
    return calim0

def test_median_XY_profile_default(create_test_data):
    calim0 = create_test_data
    
    # Call median_XY_profile function with default arguments
    result = median_XY_profile(calim0)
    
    # Validate that the result has the same shape as the input
    assert result.shape == calim0.shape, "Output shape does not match the expected shape"
    
    # Check that result contains no NaNs (since the function should fill in NaNs)
    assert not np.isnan(result).any(), "Output contains NaNs, which should be handled"

def test_median_XY_profile_no_nct(create_test_data):
    calim0 = create_test_data
    
    # Call median_XY_profile with rm_nct set to False
    result = median_XY_profile(calim0, rm_nct=False)
    
    # Validate that the result has the same shape as the input
    assert result.shape == calim0.shape, "Output shape does not match the expected shape when rm_nct is False"
    
    # Ensure the result is different when rm_nct=False compared to default (checking change of functionality)
    result_default = median_XY_profile(calim0)
    assert not np.array_equal(result, result_default), "rm_nct=False should produce a different output"

def test_median_XY_profile_margin_npixel(create_test_data):
    calim0 = create_test_data
    margin_npixel = 10
    
    # Call median_XY_profile with a different margin_npixel value
    result = median_XY_profile(calim0, margin_npixel=margin_npixel)
    
    # Validate that the result has the same shape as the input
    assert result.shape == calim0.shape, "Output shape does not match the expected shape when margin_npixel is set to 10"
    
    # Ensure the result is different compared to default margin_npixel (checking effect of margin changes)
    result_default = median_XY_profile(calim0)
    assert not np.array_equal(result, result_default), "Different margin_npixel values should produce a different output"

if __name__ == "__main__":
    pytest.main()