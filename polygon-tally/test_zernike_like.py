from zernike_like import *

import numpy as np
import pytest
import random


# testing the zernike n=0,m=0
@pytest.mark.parametrize(
    "rho, phi",
    [(0,0), (0.5, 0), (1.0, 0), (0, np.pi), (0, 2*np.pi)]
)
def test_zernike00(rho, phi):
    """asserts the zernike w/ n=0,m=0 is 1 everywhere

    rho : radial polar coordinate on the unit disk
    phi : angular polar coordinate on the unit disk
    """
    z00 = ZernikeParent(0,0)
    assert z00.zernike_nm(rho, phi, 0, 0) == 1
