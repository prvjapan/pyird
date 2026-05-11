from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
import astropy.units as u

def _get_col(table, *names):
    """Return a column by trying case variants."""
    colnames = {c.lower(): c for c in table.colnames}
    for name in names:
        key = name.lower()
        if key in colnames:
            return table[colnames[key]]
    raise KeyError(f"None of columns {names} found. Available: {table.colnames}")

def get_radec(name):
    """Get RA and Dec in degrees from an object name via SIMBAD."""

    simbad = Simbad()

    # astroquery >= 0.4.8: filter names such as 'V'
    # astroquery < 0.4.8 : 'flux(V)' style
    try:
        simbad.add_votable_fields("sp", "V", "R", "J", "H", "K")
    except Exception:
        simbad.add_votable_fields(
            "sp",
            "flux(V)", "flux(R)", "flux(J)", "flux(H)", "flux(K)"
        )

    result_table = simbad.query_object(name)
    if result_table is None or len(result_table) == 0:
        raise ValueError(f"Object not found in SIMBAD: {name}")

    ra_col = _get_col(result_table, "ra", "RA")
    dec_col = _get_col(result_table, "dec", "DEC")

    ra = ra_col[0]
    dec = dec_col[0]

    ra_unit = getattr(ra_col, "unit", None)
    dec_unit = getattr(dec_col, "unit", None)

    if ra_unit != u.deg:
        ra_unit = u.hourangle
    if dec_unit != u.deg:
        dec_unit = u.deg

    c = SkyCoord(ra, dec, unit=(ra_unit, dec_unit))
    return c.ra.degree, c.dec.degree