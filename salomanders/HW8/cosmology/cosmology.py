import numpy as np
from scipy import integrate

# speed of light
C = 299792.458  # km/s


class Cosmology(object):
    """Cosmology class implementing Cosmological Distance Functions

	import Cosmology
	cosmo = Cosmology(OmegaM=0.3, h=0.7)
	cosmo.DC(0)
	cosmo.DM(10)
	
	All equations derived from Hogg 1999 Distance Measures
	http://arxiv.org/abs/astroph/9905116
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
        return 1./np.sqrt(self.OmegaM*(1+z)**3. +
                       self.OmegaK*(1+z)**2. +
                       self.OmegaL)

    def DC(self, z):
        """Comoving Distance (Mpc)

        Computes the total line-of-sight comoving distance to a redshift of z
        """
        return self.DH*integrate.quad(self._Einv, 0, z)[0]
        
    def DM(self, z):
        """Transverse Comoving Distance (Mpc)

        The distance between two events at the 
        same redshift or distance but 
        separated on the sky by some angle dtheta
        times the transverse comoving distance
        """
        if self.OmegaK > 0:
            return self.DH/np.sqrt(self.OmegaK)*np.sinh(np.sqrt(self.OmegaK)*self.DC(z)/self.DH)
        if self.OmegaK == 0:
            return self.DC(z)
        if self.OmegaK < 0:
            return self.DH/np.sqrt(self.OmegaK)*np.sin(np.sqrt(self.OmegaK)*self.DC(z)/self.DH)
        
    def DA(self, z):
        """Angular Diameter Distance (Mpc)

        The ratio of an object's physical transverse size to 
        its angular size in radians. It is used to convert angular
        separations in telescope images to proper separations at the source.
        """
        return self.DM(z)/(1 + z)

    def DL(self, z):
        """Luminosity Distance (Mpc)

        The relationship between bolometric flux and bolometric luminosity
        """
        
        return (1 + z)*self.DM(z)

    def mu(self, z):
        """Distance Modulus (magnitudes)

        The magnitude difference between an object's observed bolometric
        flux and what it would be if it were at 10 pc
        """

        return 5.*np.log10(self.DL(z)*1e6/10.)
