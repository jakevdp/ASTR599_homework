import numpy as np
from scipy import integrate

# speed of light
C = 299792.458  # km/s


class Cosmology(object):
    """Cosmology class implementing Cosmological Distance Functions

    Parameters
    ----------
    OmegaK : float
             dimensionless curvature of space parameter
    OmegaM : float
             dimensionless matter density parameter
    OmegaL : float
             dimensionless energy density parameter
    h :  float
         dimensionless Hubble constant

    References
    ----------
    All equations taken from Hogg 00, Distance Measures in Cosmology

    """

    def __init__(self, OmegaM=0.3, h=0.7):
        # for now, we'll implement only
        # a flat cosmology for simplicity
        self.OmegaK = 0
        self.OmegaM = OmegaM
        self.OmegaL = 1. - OmegaM

        # Hubble constant, km/s/Mpc
        self.H0 = h * 100

        # Hubble Distance, Mpc
        self.DH = C / self.H0

    def _Einv(self, z):
        """implement the inverse of Eqn 14. This is the function that will be 
           integrated in order to compute other quantities
        
        Parameters
        __________
        z : float
            redshift

        Returns
        -------
        _Einv : float
                inverse of the function in Eqn 14, which is used to compute the comoving distance

        Examples
        --------
        >>> from cosmology import Cosmology
        >>> cosmo = Cosmology()
        >>> cosmo._Einv(0)
        1.0
        """

        _Einv = (1./np.sqrt(self.OmegaM * (1.+z)**3. + self.OmegaK * (1.+z)**2. + self.OmegaL))
        return _Einv

    def DC(self, z):
        """Comoving Distance (Mpc)

        Parameters
        __________
        z : float
            redshift

        Returns
        -------
        DC : float
             total line of sight comoving distance

        Examples
        --------
        >>> from cosmology import Cosmology
        >>> cosmo = Cosmology()
        >>> cosmo.DC(0)
        0.0
        """

        # Compute the comoving distance in Mpc using scipy.integrate.quad
        # following Eqn 15
        # scipy.integrate.quad returns tuple containing estimated value and upper bound - only keep estimated value
        DC = self.DH * integrate.quad(self._Einv, 0, z)[0]
        return DC

    def DM(self, z):
        """Transverse Comoving Distance (Mpc)

        Parameters
        __________
        z : float
            redshift

        Returns
        -------
        DM : float
             transverse comoving distance

        Examples
        --------
        >>> from cosmology import Cosmology
        >>> cosmo = Cosmology()
        >>> cosmo.DM(0)
        0.0
        """

        # Compute the transverse comoving distance in Mpc (Eqn 16)
        if self.OmegaK > 0. :  DM = self.DH * (1./np.sqrt(self.OmegaK)) * sinh(np.sqrt(self.OmegaK)*(self.DC(z)/self.DH))
	elif self.OmegaK == 0. : DM = self.DC(z)
	elif self.OmegaK < 0. : DM = self.DH * (1./np.sqrt(abs(self.OmegaK))) * sin(np.sqrt(abs(self.OmegaK))*(self.DC(z)/self.DH))
        return DM

    def DA(self, z):
        """Angular Diameter Distance (Mpc)

        Parameters
        __________
        z : float
            redshift

        Returns
        -------
        DA : float
             angular diameter distance

        Examples
        --------
        >>> from cosmology import Cosmology
        >>> cosmo = Cosmology()
        >>> cosmo.DA(0)
        0.0
        """

        # Compute the Angular Diameter distance in Mpc (Eqn 18)
        DA = self.DM(z) / (1.+z)
        return DA

    def DL(self, z):
        """Luminosity Distance (Mpc)

        Parameters
        __________
        z : float
            redshift

        Returns
        -------
        DL : float
             luminosity distance

        Examples
        --------
        >>> from cosmology import Cosmology
        >>> cosmo = Cosmology()
        >>> cosmo.DL(0)
        0.0
        """

        # Compute the Luminosity Distance in Mpc (Eqn 21)
        DL = (1.+z) * self.DM(z)
        return DL

    def mu(self, z):
        """Distance Modulus (magnitudes)

        Parameters
        __________
        z : float
            redshift

        Returns
        -------
        mu : float
             distance modulus

        Examples
        --------
        >>> from cosmology import Cosmology
        >>> cosmo = Cosmology()
        >>> cosmo.mu(0)
        -inf
        """

        # Compute the distance modulus (Eqn 25) -- need factor of 10**6 to account for Mpc to pc
        mu = 5. * np.log10(self.DL(z) * 10**6/ 10.)
        return mu
