import numpy as np
import pandas as pd
import pytest
from astropy.constants import c

from pyird.spec import wavrecal


def test_generate_theoretical_wavelengths_produces_expected_sequence():
    lam_low = 1000.0
    dnu = 1.0
    nline = 3

    wavelengths = wavrecal.generate_theoretical_wavelengths_of_lfc_lines(
        lam_low=lam_low,
        dnu_comb=dnu,
        nline=nline,
    )

    assert len(wavelengths) == nline + 1

    c_value = c.value
    nu0 = c_value / lam_low
    expected = [c_value / (nu0 - i * dnu) for i in range(nline + 1)]
    assert wavelengths == pytest.approx(expected)


def test_poly_lfit_recovers_coefficients_and_checks_length():
    x = np.linspace(-2, 2, 20)
    coeffs_true = np.array([0.5, -1.2, 0.7])
    y = coeffs_true[0] + coeffs_true[1] * x + coeffs_true[2] * x**2
    sig = np.ones_like(x)

    coeffs, cov, chisq, A, x_mean, x_scale = wavrecal._poly_lfit(
        x, y, sig, ma=len(coeffs_true)
    )

    fitted = wavrecal._poly_eval_scaled(x, coeffs, x_mean, x_scale)
    assert fitted == pytest.approx(y, rel=1e-6)
    assert cov.shape == (len(coeffs_true), len(coeffs_true))
    assert chisq == pytest.approx(0.0)
    assert A.shape[0] == len(x)
    assert not np.isnan(x_mean)
    assert not np.isnan(x_scale)

    short_x = np.array([0.0, 1.0])
    with pytest.raises(ValueError):
        wavrecal._poly_lfit(short_x, short_x, short_x, ma=3)


def test_reject_outliers_by_group_removes_large_residuals():
    pfit = np.arange(1.0, 11.0)
    coeffs_true = np.array([500.0, 0.2, -0.01])
    lambda_fit = coeffs_true[0] + coeffs_true[1] * pfit + coeffs_true[2] * pfit**2
    weight = np.ones_like(pfit)

    coeffs, _, _, _, x_mean, x_scale = wavrecal._poly_lfit(
        pfit, lambda_fit, weight, ma=len(coeffs_true)
    )

    lambda_with_outlier = lambda_fit.copy()
    lambda_with_outlier[8] += 5.0

    pused, lambda_used, weight_used, lambda_model, used_mask = wavrecal._reject_outliers_by_group(
        pfit,
        lambda_with_outlier,
        weight,
        x_mean,
        x_scale,
        coeffs,
        n_poly=len(coeffs_true) - 1,
        rej_lev1=3.0,
        group_size=4,
    )

    assert len(lambda_model) == len(pfit)
    assert used_mask.sum() == len(pfit) - 1
    assert not used_mask[8]
    assert len(pused) == len(lambda_used) == len(weight_used) == used_mask.sum()


def test_convert_lambda_to_pixel_interpolates_and_drops_invalid():
    lambda_obs = np.array([100.0, 101.0, 102.0, 103.0])
    lambda_peak = np.array([100.5, 102.25, 105.0, 99.5])

    pix = wavrecal._convert_lambda_to_pixel(lambda_obs, lambda_peak)

    assert len(pix) == 2
    assert pix[0] == pytest.approx(1.5)
    assert pix[1] == pytest.approx(3.25)


def test_prepare_df_fit_masks_failed_fits():
    df_obs = pd.DataFrame(
        {
            "order": [1] * 4 + [2] * 4,
            "wav": list(np.linspace(100.0, 103.0, 4)) + list(np.linspace(200.0, 203.0, 4)),
            "flux": np.arange(8, dtype=float),
        }
    )

    df_linelist = pd.DataFrame(
        {
            "order": [1, 1, 2],
            "j": [0, 1, 0],
            "x0": [100.4, 100.8, 200.5],
            "amp": [400.0, 5.0, 300.0],
            "sigma": [0.1, 0.1, 0.1],
            "fit_flg": [True, True, False],
            "lambda_fit": [100.4, 100.8, 200.5],
        }
    )

    df_fit = wavrecal._prepare_df_fit(df_obs, df_linelist)

    order1 = df_fit[df_fit["order"] == 1].reset_index(drop=True)
    assert order1.loc[0, "weight"] == pytest.approx(1 / np.sqrt(400.0))
    assert order1.loc[0, "pix"] == pytest.approx(1.4)
    assert order1.loc[1, ["pix", "weight"]].isna().all()

    order2 = df_fit[df_fit["order"] == 2].reset_index(drop=True)
    assert order2.loc[0, ["pix", "weight"]].isna().all()


def test_mk_weight_handles_negative_and_failed_lines():
    df_comblines = pd.DataFrame(
        {
            "order": [1, 1, 1, 1, 1, 1, 1],
            "j": list(range(7)),
            "x0": np.linspace(0.0, 1.0, 7),
            "amp": [4.0, 6.0, 1.0, 10.0, -5.0, 2.0, 8.0],
            "sigma": [0.1] * 7,
            "fit_flg": [True, True, False, True, True, False, True],
            "lambda_fit": np.linspace(500.0, 501.0, 7),
        }
    )

    df_weighted = wavrecal.mk_weight(df_comblines)

    assert df_weighted.loc[4, "amp"] == 0.0
    assert df_weighted.loc[2, "amp"] == pytest.approx((6.0 + 10.0) / 2)
    assert df_weighted.loc[5, "amp"] == pytest.approx(8.0)
