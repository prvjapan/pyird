import numpy as np
import pandas as pd
import pytest

from pyird.utils.irdstream import Stream1D


@pytest.fixture
def sample_stream1d(tmp_path):
    rawdir = tmp_path / "raw"
    anadir = tmp_path / "ana"
    rawdir.mkdir()
    anadir.mkdir()

    fitsid = [101, 102]
    prefix = "nw"
    extension = "mmf1"

    wav_values = np.array([100.0, 101.0])
    order_values = np.array([50, 50])
    flux_values = np.array([[10.0, 14.0], [12.0, 16.0]])
    sn_ratio_values = np.array([[50.0, 60.0], [70.0, 80.0]])
    uncertainty_values = np.array([[0.1, 0.2], [0.3, 0.4]])

    for idx, fid in enumerate(fitsid):
        df = pd.DataFrame(
            {
                "wav": wav_values,
                "order": order_values,
                "flux": flux_values[idx],
                "sn_ratio": sn_ratio_values[idx],
                "uncertainty": uncertainty_values[idx],
            }
        )
        path = rawdir / f"{prefix}{fid}{extension}.dat"
        df.to_csv(path, sep=" ", header=False, index=False)

    stream = Stream1D(
        "stream",
        rawdir,
        anadir,
        fitsid=fitsid.copy(),
        prefix=prefix,
        extension=extension,
        inst="IRD",
    )

    return stream, flux_values, uncertainty_values, wav_values


def test_specmedian_mean_computes_average_and_error(sample_stream1d):
    stream, flux_values, uncertainty_values, _ = sample_stream1d

    result = stream.specmedian(method="mean")

    expected_flux = flux_values.mean(axis=0)
    expected_err = np.sqrt((uncertainty_values**2).sum(axis=0)) / uncertainty_values.shape[0]

    assert result["flux_median"].to_numpy() == pytest.approx(expected_flux)
    assert result["flux_err"].to_numpy() == pytest.approx(expected_err)


def test_check_wavelength_range_requires_covering(sample_stream1d):
    stream, *_ = sample_stream1d

    df_ref = pd.DataFrame({"wav": [100.0, 101.0]})
    df_ok = pd.DataFrame({"wav": [99.5, 102.0]})
    stream.check_wavelength_range(df_ref, df_ok)

    df_target = pd.DataFrame({"wav": [100.2, 100.8]})
    with pytest.raises(ValueError, match="Wavelength range"):
        stream.check_wavelength_range(df_ref, df_target)


def test_recalibrate_wavelength_with_comb_validates_fiber_name(sample_stream1d):
    stream, *_ = sample_stream1d
    dummy_comb = stream.anadir / "comb.dat"

    with pytest.raises(ValueError, match="fiber must be"):
        stream.recalibrate_wavelength_with_comb(dummy_comb, fiber="science")


def test_recalibrate_wavelength_with_comb_detects_fiber_mismatch(sample_stream1d):
    stream, *_ = sample_stream1d
    dummy_comb = stream.anadir / "comb.dat"

    with pytest.raises(ValueError, match="does not match"):
        stream.recalibrate_wavelength_with_comb(dummy_comb, fiber="mmf2")
