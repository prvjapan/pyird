import astropy.units as u
import pytest
from astropy.table import Table

import pyird.utils.getinfo as getinfo


class FakeSimbad:
    response = None
    added_fields = []
    queried_names = []

    def add_votable_fields(self, *fields):
        self.added_fields.append(fields)

    def query_object(self, name):
        self.queried_names.append(name)
        return self.response


@pytest.fixture(autouse=True)
def fake_simbad(monkeypatch):
    FakeSimbad.response = None
    FakeSimbad.added_fields = []
    FakeSimbad.queried_names = []
    monkeypatch.setattr(getinfo, "Simbad", FakeSimbad)


def test_get_radec_from_degree_columns():
    table = Table({"ra": [346.86964875], "dec": [21.134252777777778]})
    table["ra"].unit = u.deg
    table["dec"].unit = u.deg
    FakeSimbad.response = table

    ra, dec = getinfo.get_radec("HR8799 A")

    assert ra == pytest.approx(346.86964875)
    assert dec == pytest.approx(21.134252777777778)
    assert FakeSimbad.queried_names == ["HR8799 A"]
    assert FakeSimbad.added_fields == [("sp", "V", "R", "J", "H", "K")]


def test_get_radec_from_sexagesimal_columns():
    FakeSimbad.response = Table(
        {"RA": ["23 07 28.7157"], "DEC": ["+21 08 03.3100"]}
    )

    ra, dec = getinfo.get_radec("HR8799 A")

    assert ra == pytest.approx(346.86964875, abs=1e-8)
    assert dec == pytest.approx(21.13425278, abs=1e-8)


@pytest.mark.parametrize("response", [None, Table()])
def test_get_radec_raises_for_missing_object(response):
    FakeSimbad.response = response

    with pytest.raises(ValueError, match="Object not found in SIMBAD: unknown"):
        getinfo.get_radec("unknown")


def test_get_radec_uses_legacy_field_names_when_needed(monkeypatch):
    class LegacySimbad(FakeSimbad):
        def add_votable_fields(self, *fields):
            self.added_fields.append(fields)
            if "V" in fields:
                raise ValueError("modern field names are not supported")

    table = Table({"ra": [1.0], "dec": [2.0]})
    table["ra"].unit = u.deg
    table["dec"].unit = u.deg
    LegacySimbad.response = table
    monkeypatch.setattr(getinfo, "Simbad", LegacySimbad)

    getinfo.get_radec("target")

    assert LegacySimbad.added_fields == [
        ("sp", "V", "R", "J", "H", "K"),
        ("sp", "flux(V)", "flux(R)", "flux(J)", "flux(H)", "flux(K)"),
    ]
