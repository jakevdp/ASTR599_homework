<<<<<<< HEAD
#%%file cosmology.py
=======
>>>>>>> 03749cc8c76066974d9fc99a7bf0963e4e82b7f8
import numpy as np
from scipy import integrate

# speed of light
C = 299792.458  # km/s


class Cosmology(object):
    """Cosmology class implementing Cosmological Distance Functions

    Equation used in this code are taken from Hogg 1999
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
<<<<<<< HEAD
	return 1./(self.OmegaM*(1+z)**3.+self.OmegaK*(1+z)**2.+self.OmegaL)**(0.5)      
=======
	return 1./(self.OmegaM*(1+z)**3.+self.OmegaK*(1+z)**2.+self.OmegaL)**(0.5)       
>>>>>>> 03749cc8c76066974d9fc99a7bf0963e4e82b7f8
	

    def DC(self, z):
        """Comoving Distance (Mpc)

        Total line of sight comoving distance
	Input z        
	"""
        # Compute the comoving distance in Mpc using scipy.integrate.quad
        # following Eqn 15
	from scipy import integrate
	x,_ = integrate.quad(self._Einv,0,z)
	return self.DH*x

    def DM(self, z):
        """Transverse Comoving Distance (Mpc)

        Comoving distance between two events at the same redshit but separated on the sky by some angle dtheta is DMdtheta
	Input z
        """
        # Compute the transverse comoving distance in Mpc (Eqn 16)
	return self.DC(z)

    def DA(self, z):
        """Angular Diameter Distance (Mpc)

        Ratio of an object's physical transverse size to its angular size
	Input z
        """
        # Compute the Angular Diameter distance in Mpc (Eqn 18)
	return self.DM(z)/(1.+z)
	
    def DL(self, z):
        """Luminosity Distance (Mpc)

        Relationship between bolometric flux and bolometric luminosity
	Input z
        """
        # Compute the Luminosity Distance in Mpc (Eqn 21)
	return (1.+z)**2.*self.DA(z)
	

    def mu(self, z):
        """Distance Modulus (magnitudes)

        Input z
        """
        # Compute the distance modulus (Eqn 25)
<<<<<<< HEAD
	return 5.*np.log10(1000000.*self.DL(z)/10.)
=======
	return 5.*numpy.log10(self.DL(z)/10.)
>>>>>>> 03749cc8c76066974d9fc99a7bf0963e4e82b7f8
