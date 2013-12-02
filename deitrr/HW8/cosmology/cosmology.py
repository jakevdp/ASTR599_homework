import numpy as np
from scipy import integrate

# speed of light
C = 299792.458  # km/s


class Cosmology(object):
    """Cosmology class implementing Cosmological Distance Functions

    Creates a class object with cosmological distances and parameters
    as its methods.  Based on formulae of David Hogg (1999).
    
    Parameters
    ----------
    OmegaM : float, optional
        Matter density parameter. Default = 0.3
    h : float, optional
        dimensionless Hubble constant. Default = 0.7
        
    References
    ----------
    Hogg, D., 1999. arXiv:astro-ph/9905116v4
    """
    def __init__(self, OmegaM=0.3, h=0.7):
        # for now, we'll implement only
        # a flat cosmology for simplicity
        self.OmegaK = 0.0
        self.OmegaM = OmegaM
        self.OmegaL = 1. - OmegaM

        # Hubble constant, km/s/Mpc
        self.H0 = h * 100

        # Hubble Distance, Mpc
        self.DH = C / self.H0

    def _Einv(self, z):
        """Returns the inverse of Hogg equation 14 
        
        Parameters
        ----------
        z : float
            redshift
        
        Returns
        -------
        y : float
            E^(-1) where E is given by Hogg eqn 14
        """
        return 1.0/np.sqrt(self.OmegaM*(1.0+z)**3.0 + self.OmegaK \
                        *(1.0+z)**2.0 + self.OmegaL)

    def DC(self, z):
        """Comoving Distance (Mpc)

        Parameters
        ----------
        z : float
            redshift
        
        Returns
        -------
        y : float
            The comoving distance in Mpc, given by Hogg eqn 15
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DC(1.0)
        3303.8288058874678
        """
        # Compute the comoving distance in Mpc using scipy.integrate.quad
        # following Eqn 15
        return self.DH * integrate.quad(self._Einv,0,z)[0]

    def DM(self, z):
        """Transverse Comoving Distance (Mpc)

        Parameters
        ----------
        z : float
            redshift
        
        Returns
        -------
        y : float
            The transverse comoving distance in Mpc, given by Hogg eqn 16
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DM(1.0)
        3303.8288058874678
        """
        # Compute the transverse comoving distance in Mpc (Eqn 16)
        if self.OmegaK > 0.0:
            return self.DH / np.sqrt(self.OmegaK) * \
                    np.sinh(np.sqrt(self.OmegaK)*self.DC(z)/self.DH)
        elif self.OmegaK == 0.0:
            return self.DC(z)
        elif self.OmegaK < 0.0:
            return self.DH / np.sqrt(np.abs(self.OmegaK)) * \
                    np.sin(np.sqrt(np.abs(self.OmegaK))*self.DC(z)/self.DH)

    def DA(self, z):
        """Angular Diameter Distance (Mpc)
        
        Parameters
        ----------
        z : float
            redshift
        
        Returns
        -------
        y : float
            The angular diameter distance in Mpc, given by Hogg eqn 18
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DA(1.0)
        1651.9144029437339
        """
        # Compute the Angular Diameter distance in Mpc (Eqn 18)
        return self.DM(z) / (1.0+z)

    def DL(self, z):
        """Luminosity Distance (Mpc)
        
        Parameters
        ----------
        z : float
            redshift
        
        Returns
        -------
        y : float
            The luminosity distance in Mpc, given by Hogg eqn 21
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DL(1.0)
        6607.6576117749355
        """
        # Compute the Luminosity Distance in Mpc (Eqn 21)
        return (1.0+z)*self.DM(z)

    def mu(self, z):
        """Distance Modulus (magnitudes)

        Parameters
        ----------
        z : float
            redshift
        
        Returns
        -------
        y : float
            The distance modulus, given by Hogg eqn 25
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.mu(1.0)
        44.100237655543722
        """
        # Compute the distance modulus (Eqn 25)
        return 5.0*np.log10(self.DL(z)*(1.e6)/10.0)