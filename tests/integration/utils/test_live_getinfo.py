import os

import pytest
from pyvo.dal.exceptions import DALServiceError
from requests.exceptions import ConnectionError, Timeout

from pyird.utils.getinfo import get_radec


pytestmark = [
    pytest.mark.online,
    pytest.mark.simbad,
    pytest.mark.skipif(
        os.environ.get("PYIRD_RUN_ONLINE_TESTS") != "1",
        reason="Set PYIRD_RUN_ONLINE_TESTS=1 to run live SIMBAD tests",
    ),
]


def test_get_radec_with_live_simbad():
    try:
        ra, dec = get_radec("HR8799 A")
    except (ConnectionError, Timeout, DALServiceError) as error:
        pytest.xfail(f"SIMBAD is unavailable: {error}")

    assert ra == pytest.approx(346.86964875, abs=1e-5)
    assert dec == pytest.approx(21.13425278, abs=1e-5)
