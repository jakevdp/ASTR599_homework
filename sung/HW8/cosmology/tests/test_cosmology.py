from numpy.testing import assert_allclose
from .. import Cosmology


# [OmegaM, h, z, DC]
DC_RESULTS =\
[[0.2, 0.6, 0, 0.0],
 [0.2, 0.6, 0.5, 2285.61],
 [0.2, 0.6, 1.0, 4114.6],
 [0.2, 0.7, 0, 0.0],
 [0.2, 0.7, 0.5, 1959.09],
 [0.2, 0.7, 1.0, 3526.8],
 [0.3, 0.6, 0, 0.0],
 [0.3, 0.6, 0.5, 2203.4],
 [0.3, 0.6, 1.0, 3854.47],
 [0.3, 0.7, 0, 0.0],
 [0.3, 0.7, 0.5, 1888.63],
 [0.3, 0.7, 1.0, 3303.83],
 [0.4, 0.6, 0, 0.0],
 [0.4, 0.6, 0.5, 2131.92],
 [0.4, 0.6, 1.0, 3649.21],
 [0.4, 0.7, 0, 0.0],
 [0.4, 0.7, 0.5, 1827.36],
 [0.4, 0.7, 1.0, 3127.89]]


# [OmegaM, h, z, DM]
DM_RESULTS =\
[[0.2, 0.6, 0, 0.0],
 [0.2, 0.6, 0.5, 2285.61],
 [0.2, 0.6, 1.0, 4114.6],
 [0.2, 0.7, 0, 0.0],
 [0.2, 0.7, 0.5, 1959.09],
 [0.2, 0.7, 1.0, 3526.8],
 [0.3, 0.6, 0, 0.0],
 [0.3, 0.6, 0.5, 2203.4],
 [0.3, 0.6, 1.0, 3854.47],
 [0.3, 0.7, 0, 0.0],
 [0.3, 0.7, 0.5, 1888.63],
 [0.3, 0.7, 1.0, 3303.83],
 [0.4, 0.6, 0, 0.0],
 [0.4, 0.6, 0.5, 2131.92],
 [0.4, 0.6, 1.0, 3649.21],
 [0.4, 0.7, 0, 0.0],
 [0.4, 0.7, 0.5, 1827.36],
 [0.4, 0.7, 1.0, 3127.89]]


# [OmegaM, h, z, DA]
DA_RESULTS =\
[[0.2, 0.6, 0, 0.0],
 [0.2, 0.6, 0.5, 1523.74],
 [0.2, 0.6, 1.0, 2057.3],
 [0.2, 0.7, 0, 0.0],
 [0.2, 0.7, 0.5, 1306.06],
 [0.2, 0.7, 1.0, 1763.4],
 [0.3, 0.6, 0, 0.0],
 [0.3, 0.6, 0.5, 1468.93],
 [0.3, 0.6, 1.0, 1927.23],
 [0.3, 0.7, 0, 0.0],
 [0.3, 0.7, 0.5, 1259.08],
 [0.3, 0.7, 1.0, 1651.91],
 [0.4, 0.6, 0, 0.0],
 [0.4, 0.6, 0.5, 1421.28],
 [0.4, 0.6, 1.0, 1824.61],
 [0.4, 0.7, 0, 0.0],
 [0.4, 0.7, 0.5, 1218.24],
 [0.4, 0.7, 1.0, 1563.95]]

# [OmegaM, h, z, DL]
DL_RESULTS =\
[[0.2, 0.6, 0, 0.0],
 [0.2, 0.6, 0.5, 3428.41],
 [0.2, 0.6, 1.0, 8229.21],
 [0.2, 0.7, 0, 0.0],
 [0.2, 0.7, 0.5, 2938.64],
 [0.2, 0.7, 1.0, 7053.61],
 [0.3, 0.6, 0, 0.0],
 [0.3, 0.6, 0.5, 3305.09],
 [0.3, 0.6, 1.0, 7708.93],
 [0.3, 0.7, 0, 0.0],
 [0.3, 0.7, 0.5, 2832.94],
 [0.3, 0.7, 1.0, 6607.66],
 [0.4, 0.6, 0, 0.0],
 [0.4, 0.6, 0.5, 3197.88],
 [0.4, 0.6, 1.0, 7298.42],
 [0.4, 0.7, 0, 0.0],
 [0.4, 0.7, 0.5, 2741.04],
 [0.4, 0.7, 1.0, 6255.79]]

# [OmegaM, h, z, mu]
MU_RESULTS = \
[[0.2, 0.6, 0.1, 38.666337868889762],
 [0.2, 0.6, 0.5, 42.675462188569789],
 [0.2, 0.6, 1.0, 44.57678993423005],
 [0.2, 0.7, 0.1, 38.331603920736697],
 [0.2, 0.7, 0.5, 42.340728240416723],
 [0.2, 0.7, 1.0, 44.242055986076984],
 [0.3, 0.6, 0.1, 38.649938522544012],
 [0.3, 0.6, 0.5, 42.595919369693959],
 [0.3, 0.6, 1.0, 44.434971603696795],
 [0.3, 0.7, 0.1, 38.31520457439094],
 [0.3, 0.7, 0.5, 42.261185421540894],
 [0.3, 0.7, 1.0, 44.100237655543722],
 [0.4, 0.6, 0.1, 38.633904578021593],
 [0.4, 0.6, 0.5, 42.524308199311712],
 [0.4, 0.6, 1.0, 44.316144393006496],
 [0.4, 0.7, 0.1, 38.299170629868527],
 [0.4, 0.7, 0.5, 42.189574251158646],
 [0.4, 0.7, 1.0, 43.981410444853431]]


def assert_equal_to_2_decimals(x, y):
    assert_allclose(x, y, atol=0.01)


def test_DC():
    """Test Comoving Distance"""
    for (OmegaM, h, z, DC) in DC_RESULTS:
        cosmo = Cosmology(OmegaM=OmegaM, h=h)
        yield assert_equal_to_2_decimals, DC, cosmo.DC(z)


def test_DM():
    """Test Transverse Comoving Distance"""
    for (OmegaM, h, z, DM) in DM_RESULTS:
        cosmo = Cosmology(OmegaM=OmegaM, h=h)
        yield assert_equal_to_2_decimals, DM, cosmo.DM(z)


def test_DA():
    """Test Angular Diameter Distance"""
    for (OmegaM, h, z, DA) in DA_RESULTS:
        cosmo = Cosmology(OmegaM=OmegaM, h=h)
        yield assert_equal_to_2_decimals, DA, cosmo.DA(z)


def test_DL():
    """Test Luminosity Distance"""
    for (OmegaM, h, z, DL) in DL_RESULTS:
        cosmo = Cosmology(OmegaM=OmegaM, h=h)
        yield assert_equal_to_2_decimals, DL, cosmo.DL(z)


def test_mu():
    """Test Distance Modulus"""
    for (OmegaM, h, z, MU) in MU_RESULTS:
        cosmo = Cosmology(OmegaM=OmegaM, h=h)
        yield assert_equal_to_2_decimals, MU, cosmo.mu(z)
