import numpy as np
from scipy import integrate

# speed of light
C = 299792.458  # km/s


class Cosmology(object):
    """

    Cosmology class implementing Cosmological Distance Functions

    Methods
    -------
    DC - computes the comoving distance
    DM - computes the transverse comoving distance
    DA - computes the angular diameter distance
    DL - computes the luminosity distance
    mu - computes the distance modulus
    
    Notes
    -----
    See methods for details for calculating different cosmological distances
    
    References
    ----------
    Higgs (1999)
    
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
        # implement the inverse of Eqn 14. This is the function that will be
        # integrated in order to compute other quantities
        return 1.0 / np.sqrt(self.OmegaM*((1.+z)**3) + self.OmegaK*((1.+z)**2) + self.OmegaL)

    def DC(self, z):
        """
        Computes the Comoving Distance (Mpc)

        Parameters
        ----------
        z : float 
            The redshift z which sets the upper limit of the integral
        
        Returns
        -------
        y : float
            The corresponding co-moving distance in Mpc
            
        Examples
        --------
        >>> cosmo=Cosmology()
        >>> cosmo.DC(1.2)
        3763.4313597326077
        
        Notes
        -----
        The integral is evaluated using scipy.integrate.quad
        
        References
        ----------
        Eqn. 15 of Higgs (1999)
        """
        # Compute the comoving distance in Mpc using scipy.integrate.quad
        # following Eqn 15
        #Ez = lambda z: (self.OmegaM*(1.+z)**3 + self.OmegaK*(1.+z)**2 + self.OmegaL)**(-0.5)
        
        y = self.DH * integrate.quad(self._Einv, 0., float(z))[0]
        return y

    def DM(self, z):
        """
        Computes the Transverse Comoving Distance (Mpc)

        Parameters
        ----------
        z : float 
            The redshift z at which the distance is to be evaluated
        
        Returns
        -------
        y : float
            The corresponding transverse co-moving distance in Mpc
            
        Examples
        --------
        >>> cosmo=Cosmology()
        >>> cosmo.DM(1.2)
        3763.4313597326077
        
        Notes
        -----
        The returned value is different for omegaK > 0, omegaK = 0, and omegaK < 0.
        
        References
        ----------
        Eqn. 16 of Higgs (1999)
        """
        # Compute the transverse comoving distance in Mpc (Eqn 16)
        if self.OmegaK > 0:
            y = self.DH * (self.OmegaK)**(-0.5) * np.sinh(self.OmegaK**0.5 * self.DC(z)/self.DH)
        if self.OmegaK == 0:
            y = self.DC(z)
        if self.OmegaK < 0:
            y = self.DH * abs(self.OmegaK)**(-0.5) * np.sinh(abs(self.OmegaK)**0.5 * self.DC(z) / self.DH)
        return y

    def DA(self, z):
        """
        Computes the Angular Diameter Distance (Mpc)

        Parameters
        ----------
        z : float 
            The redshift z at which the distance is to be evaluated
        
        Returns
        -------
        y : float
            The corresponding angular diameter distance in Mpc
            
        Examples
        --------
        >>> cosmo=Cosmology()
        >>> cosmo.DA(1.2)
        1710.650618060276
        
        Notes
        -----

        References
        ----------
        Eqn. 18 of Higgs (1999)
        """
        # Compute the Angular Diameter distance in Mpc (Eqn 18)
        y = self.DM(z) / (1. + z)
        return y
        
    def DL(self, z):
        """
        Computes the Luminosity Distance (Mpc)

        Parameters
        ----------
        z : float 
            The redshift z at which the distance is to be evaluated
        
        Returns
        -------
        y : float
            The corresponding luminosity distance in Mpc
            
        Examples
        --------
        >>> cosmo=Cosmology()
        >>> cosmo.DL(1.2)
        8279.548991411737
        
        Notes
        -----
        DL may be evaluated either (1+z) * DM or (1+z)^2 * DA. This uses the DM.

        References
        ----------
        Eqn. 21 of Higgs (1999)
        """
        # Compute the Luminosity Distance in Mpc (Eqn 21)
        y = (1. + z) * self.DM(z)
        return y

    def mu(self, z):
        """
        Computes the Distance Modulus (magnitudes)

        Parameters
        ----------
        z : float 
            The redshift z at which the distance modulus is evaluated
        
        Returns
        -------
        y : float
            The corresponding distance modulus in magnitudes
            
        Examples
        --------
        >>> cosmo=Cosmology()
        >>> cosmo.mu(1.2)
        44.590033401390663
        
        Notes
        -----
        The returned value is different for omegaK > 0, omegaK = 0, and omegaK < 0.
        
        References
        ----------
        Eqn. 25 of Higgs (1999)
        """
        # Compute the distance modulus (Eqn 25)
        y = 5.0 * np.log10(self.DL(z) * 1e6 / 10.)
        return y
