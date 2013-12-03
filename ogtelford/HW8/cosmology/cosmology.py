import numpy as np
from scipy import integrate

# speed of light
C = 299792.458  # km/s


class Cosmology(object):
    """Cosmology class implementing Cosmological Distance Functions
    
    The functions within this class calculate the comoving distance,
    transverse comoving distance, angular diameter distance, luminosity distance,
    and distance modulus at a given redshift.
    
    Parameters
    ----------
    OmegaM : float, optional
             The dimensionless mass density parameter. 0.3 by default.
    
    h : float, optional
        A dimensionless quantity parameterizing the Hubble constant. 0.7 by default.
    

    References
    ----------
    ..[1] Hogg, David W., "Distance Measures in Cosmology," arXiv, 1999.
    
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
        
        return 1./np.sqrt(self.OmegaM*((1+z)**3) + self.OmegaK*((1+z)**2) + self.OmegaL)

    def DC(self, z):
        """Comoving Distance (Mpc)
        
        Returns the comoving distance in units of Mpc at a given redshift.

        Parameters
        ----------
        z : float
            the cosmological redshift. 
    
        Returns
        -------
        DC : float
            the comoving distance in Mpc at the redshift z as defined in eqn 15 (Hogg).
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DC(1)
        3303.8288058874678
        """
        # Compute the comoving distance in Mpc using scipy.integrate.quad
        # following Eqn 15
        
        return self.DH*integrate.quad(self._Einv, 0, z)[0] # is this self._Einv the correct syntax?

    def DM(self, z):
        """Transverse Comoving Distance (Mpc)
        
        Returns the transverse comoving distance in units of Mpc at a given redshift.

        Parameters
        ----------
        z : float
            the cosmological redshift. 
    
        Returns
        -------
        DM : float
            the transverse comoving distance in Mpc at the redshift z as defined in eqn 16 (Hogg).
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DM(1)  
        3303.8288058874678
        """
        # Compute the transverse comoving distance in Mpc (Eqn 16)
        
        if self.OmegaK == 0:
            DM = self.DC(z)
        else:
            DM = (self.DH/np.sqrt(np.abs(self.OmegaK)))*np.sinh((np.sqrt(np.abs(self.OmegaK))/self.DH)*self.DC(z))
        return DM
            
    def DA(self, z):
        """Angular Diameter Distance (Mpc)

        Returns the angular diameter distance in units of Mpc at a given redshift.

        Parameters
        ----------
        z : float
            the cosmological redshift. 
    
        Returns
        -------
        DA : float
            the angular diameter distance in Mpc at the redshift z as defined in eqn 18 (Hogg).
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DA(1)
        1651.9144029437339
        """
        # Compute the Angular Diameter distance in Mpc (Eqn 18)
        
        return self.DM(z)/(1+z)

    def DL(self, z):
        """Luminosity Distance (Mpc)

        Returns the luminosity distance in units of Mpc at a given redshift.

        Parameters
        ----------
        z : float
            the cosmological redshift. 
    
        Returns
        -------
        DL : float
            the luminosity distance in Mpc at the redshift z as defined in eqn 21 (Hogg).
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.DL(1)
        6607.6576117749355
        """
        # Compute the Luminosity Distance in Mpc (Eqn 21)
        
        return self.DM(z)*(1+z)

    def mu(self, z):
        """Distance Modulus (magnitudes)
        
        Returns the distance modulus in units of magnitudes at a given redshift.
        The distance modulus is the difference between the object's observed
        bolometric magnitude and the bolometric magnitude it would have if it
        were at a distance of 10 pc.

        Parameters
        ----------
        z : float
            the cosmological redshift. 
    
        Returns
        -------
        mu : float
            the distance modulus in magnitudes at the redshift z as defined in eqn 25 (Hogg).
            
        Examples
        --------
        >>> cosmo = Cosmology()
        >>> cosmo.mu(1)
        44.100237655543722
        """
        # Compute the distance modulus (Eqn 25)
        
        return 5*np.log10(self.DL(z)*(1e6)/10)