import pytest
import pandas as pd
import numpy as np
from pyird.cleaning.remove_airglow import load_oliva15, wav_around_airglow, df_mask_airglow

@pytest.fixture
def test_data():
    # Prepare sample data for testing
    wav_obs = pd.Series([1500, 1501, 1502, 1503, 1504, 1505, 1506])
    flux = pd.Series([0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07])
    df_obs = pd.DataFrame({"wav": wav_obs, "flux": flux})
    return df_obs

def test_load_oliva15():
    # Test if load_oliva15 loads the wavelengths correctly and sorts them
    wav_all = load_oliva15()
    
    # Check the type and content of the result
    assert isinstance(wav_all, pd.Series)
    assert wav_all.is_monotonic_increasing  # Should be sorted in increasing order
    assert np.all(wav_all >= 0)  # Ensure that no negative wavelengths are present

def test_wav_around_airglow():
    # Test the wav_around_airglow function
    wav_obs = pd.Series([1500, 1501, 1502, 1503, 1504, 1505, 1506])
    wav_airglow = [1502, 1505]
    ind_mask = wav_around_airglow(wav_obs, wav_airglow, mask_width=0.035)
    
    # Check if the indices are correct
    assert len(ind_mask) == 2  # Two airglow features
    
    # Verify the indices that are masked (1502 and 1505 should mask surrounding points)
    assert 2 in ind_mask  # 1502 should mask the point at index 2
    assert 5 in ind_mask  # 1505 should mask the point at index 5

def test_df_mask_airglow(test_data):
    # Test the df_mask_airglow function
    df_obs = test_data
    wav_airglow = [1502, 1505]
    ind_mask = wav_around_airglow(df_obs['wav'], wav_airglow, mask_width=0.035)
    
    # Mask the airglow
    df_masked = df_mask_airglow(df_obs, ind_mask, plot=False)
    
    # Check if the number of rows in df_masked is reduced after masking
    assert len(df_masked) == len(df_obs) - len(ind_mask)
    
    # Ensure that the masked indices are removed
    assert not any(df_masked['wav'].isin(df_obs['wav'][ind_mask]))  # Should not contain masked wavelengths

@pytest.mark.parametrize("wav_airglow, expected_mask_count", [
    ([1502], 1),  # Test with a single airglow feature
    ([1502, 1505], 2)  # Test with two airglow features
])
def test_wav_around_airglow_parametrized(wav_airglow, expected_mask_count):
    # Parametrized test case for different airglow positions
    wav_obs = pd.Series([1500, 1501, 1502, 1503, 1504, 1505, 1506])
    ind_mask = wav_around_airglow(wav_obs, wav_airglow, mask_width=0.035)
    
    # Check if the number of indices masked matches the expected count
    assert len(ind_mask) == expected_mask_count
