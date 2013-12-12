import numpy as np
from scipy import integrate

# speed of light
C = 299792.458  # km/s


class Cosmology(object):
    """Cosmology class implementing Cosmological Distance Functions

    This module implements a class to calculate cosmological paraments given
    certain observables. The equations are implemented for a flat cosmology only.
    
    Parameters
    ----------
    OmegaM : float
        Mass density parameter
        
    h : float
        Hubble parameter where the Hubble constant H0=h*100. km/s/Mpc
    
    H0 : float
        Hubble constant. See above.
    
    OmegaK : float
        curvature density parameter
        
    OmegaL : float
        Lambda density parameter (dark energy/cosmological constant)
    
    DH : float
        Hubble Distance in Mpc
    
    References
    ----------
    The equations are obtrained from Hogg, David W. 2000, arXiv:astro-ph/9905116
    
    """
    def __init__(self, OmegaM=0.3, h=0.7):
        # for now, we'll implement only
        # a flat cosmology for simplicity
        self.OmegaK = 0.
        self.OmegaM = OmegaM
        self.OmegaL = 1. - OmegaM

        # Hubble constant, km/s/Mpc
        self.H0 = h * 100.

        # Hubble Distance, Mpc
        self.DH = C / self.H0

    def Einv(self, z):
        # implement the inverse of Eqn 14. This is the function that will be
        # integrated in order to compute other quantities
        return 1./(np.sqrt(self.OmegaM*(1.0+z)*(1.0+z)*(1.0+z) + self.OmegaK*(1.0+z)*(1.0+z) + self.OmegaL))

    def DC(self, z):
        """
        Compute the Comoving Distance in Mpc
    
        Parameters
        ----------
        z : float
        Redshift
    
        Returns
        -------
        self.DC : float
            DC = DH*integral(Einv,(0,z))
        
        Notes
        -------
        From equation 15 in Hogg 2000
        
        Examples
        --------
        >>> Cosmology().DC(0.0)
        0.0
        >>> Cosmology().DC(0.1)
        418.4544876277074
        
        """
        # Compute the comoving distance in Mpc using scipy.integrate.quad
        # following Eqn 15 
        return self.DH*integrate.quad(self.Einv,0.0,z)[0]
    
    def DM(self, z):
        """
        Compute the Transverse Comoving Distance in Mpc

        Parameters
        ----------
        z : float
        Redshift
    
        Returns
        -------
        self.DM : float
            DM = DC(z)
        
        Notes
        -------
        This is trivial since we are implementing in a flat cosmology.
        In this case DM(z)=DC(Z). From equation 16 in Hogg 2000.
        
        Examples
        --------
        >>> Cosmology().DM(0.0)
        0.0
        >>> Cosmology().DM(0.1)
        418.4544876277074
    
        """
        # Compute the transverse comoving distance in Mpc (Eqn 16)
        #In a flat cosmology this is equal to the comoving distance
        return self.DC(z)

    def DA(self, z):
        """
        Compute the Angular Diameter Distance in Mpc

        Parameters
        ----------
        z : float
        Redshift
    
        Returns
        -------
        self.DA : float
            DA = DM(z)/(1+z)
        
        Notes
        -------
        The angular diameter distance is the ratio of an object's
        physical size to its angular size on the sky in radians.
        From equation 18 in Hogg 2000.
        
        Examples
        --------
        >>> Cosmology().DA(0.0)
        0.0
        >>> Cosmology().DA(0.1)
        380.41317057064305
    
        """
        # Compute the Angular Diameter distance in Mpc (Eqn 18)
        return self.DM(z)/(1+z)
        
    def DL(self, z):
        """
        Computes the Luminosity Distance in Mpc

        Parameters
        ----------
        z : float
        Redshift
    
        Returns
        -------
        self.DL : float
            DL = DM(z) * (1+z)
        
        Notes
        -------
        The luminosity distance is related to the observed
        integrated flux from an object compared to its intrinsic luminosity.
        For objects at low redshifts this is equal to the 'true' spatial distance.
        From equation 21 in Hogg 2000.
        
        Examples
        --------
        >>> Cosmology().DL(0.0)
        0.0
        >>> Cosmology().DL(0.1)
        460.2999363904782
        
        """
        # Compute the Luminosity Distance in Mpc (Eqn 21)
        return self.DM(z)*(1+z)
        
    def mu(self, z):
        """
        Computes the Distance Modulus in magnitudes

        Parameters
        ----------
        z : float
        Redshift
    
        Returns
        -------
        self.mu : float
            mu = 5. * log10(DL(z)/10.) + 5*log10(h)
        
        Notes
        -------
        The distance modulus is used to relate the apparent
        magnitude to the absolute magnitude of the object.
        From equation 25 in Hogg 2000.
        
        Examples
        --------
        >>> Cosmology().mu(0.1)
        38.31520457439094
        """
        # Compute the distance modulus (Eqn 25)
        return 5. * np.log10(self.DL(z)/0.00001) 
        #bleh, converted from pc to mpc above