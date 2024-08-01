import numpy as np
from math import factorial

class ZernikeParent:
    """Parent class for the ZBasis, KBasis, HBasis classes""" 

    # initialization
    def __init__(self, n, m):
        """
        n : Zernike order
        m : Zernike sub-order
        """
        
        # checking validity
        self.n_checker(n)
        self.m_checker(n, m)
        
        # initalizing parameters
        self.n = n
        self.m = m


    # initialization checkers 
    def n_checker(self, n):
        """Verfies the n associated with the Zernike order is valid.

        n : Zernike order
        """
        if type(n) is not int:
            raise TypeError("invalid n; make sure n is an integer")

    def m_checker(self, n, m):
        """Verifies the m associated with the Zernike order is valid for the given n.
        
        n : Zernike order
        m : Zernike sub-order
        """
        if m not in np.arange(-n, n+1, 2):
            raise ValueError("invalid m; enter an m from [-n, -n+2, -n+4, ... , n]")

    
    # zernike intermediate functions
    def n_nm(self, n, m):
        """Calculates the normalization constant for the Zernike order

        n : Zernike order
        m : Zernike sub-order
        """
        if m != 0:
            return (2 * (n + 1))**(1 / 2)
        else:
            return (n + 1)**(1 / 2)
    

    def r_nm(self, rho, n, m):
        """Calculates an intermediate function needed for the Zernike polynomials
        
        rho : radial polar cordinate on the unit disk
        n : Zernike order
        m : Zernike sub-order
        """

        # variables
        max_index = int((n - abs(m)) / 2)
        R_nm = 0

        # summation
        for k in range(max_index + 1):
            numerator = (-1)**k * factorial(n - k) * rho**(n - 2 * k)
            denominator = factorial(k) * factorial(int((n + abs(m)) / 2 - k)) * factorial(int((n - abs(m)) / 2 - k))

            R_nm += numerator / denominator

        return R_nm


    def zernike_nm(self, rho, phi, n, m):
        """creates the zernike function associated with the respective n and m values

        rho : radial polar coordinate on the unit disk
        phi : angular polar coordinate on the unit disk
        n : Zernike order
        m : Zernike sub-order
        """
        # checking
        if rho > 1.0 or rho < 0.0:
            raise ValueError("rho out of bounds; rho ranges from [0,1]")
        
        if phi > 2 * np.pi or phi < 0.0:
            raise ValueError("phi out of bounds; phi ranges from [0,2pi)")


        # calculating the value of the Zernike polynomial at a certain coordinate
        if m >= 0:
            return self.n_nm(n, m) * self.r_nm(rho, n, m) * np.cos(m * phi)
        else:
            return -1 * self.n_nm(n, m) * self.r_nm(rho, n, m) * np.sin(m * phi)


class ZBasis(ZernikeParent):
    """class for the zernike basis vectors, which exist on the unit disk"""





# outputs
z1 = ZernikeParent(0,0)
print(z1.r_nm(.5,0,0))
print(z1.zernike_nm(0,0,0,0))
