"""get information."""

from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy import units as u
from packaging.version import parse
import astroquery


def get_radec(name):
    """Get ra and dec from name via Simbad.

    Args:
        name: object name

    Returns:
        ra in degree
        dec in degree
    """

    Simbad.SIMBAD_URL = 'http://simbad.u-strasbg.fr/simbad/sim-script'
    if parse(astroquery.__version__) <= parse("0.4.7"):
        Simbad.add_votable_fields(
            'sp', 'flux(V)', 'flux(R)', 'flux(J)', 'flux(H)', 'flux(K)')
    else:
        Simbad.add_votable_fields('sp', 'flux')

    result_table = Simbad.query_object(name)
    namex = result_table['MAIN_ID'][0]  # .decode('utf-8')
    ra = result_table['RA'][0]
    dec = result_table['DEC'][0]
    c = SkyCoord(str(ra)+' '+str(dec), unit=(u.hourangle, u.deg))
    ra = c.ra.degree
    dec = c.dec.degree
    return ra, dec
