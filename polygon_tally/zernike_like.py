import numpy as np
from math import factorial

class ZernikeParent:
    """Parent class for the ZBasis, KBasis, HBasis classes""" 

    # initialization
    def __init__(self, n, m):
        """
        Parameters
        ----------
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

        Parameters
        ----------
        n : Zernike order
        """
        if type(n) is not int:
            raise TypeError("invalid n; make sure n is an integer")

    def m_checker(self, n, m):
        """Verifies the m associated with the Zernike order is valid for the given n.
        
        Parameters
        ----------
        n : Zernike order
        m : Zernike sub-order
        """
        if m not in np.arange(-n, n+1, 2):
            raise ValueError("invalid m; enter an m from [-n, -n+2, -n+4, ... , n]")

    
    # zernike intermediate functions
    def n_nm(self, n, m):
        """Calculates the normalization constant for the Zernike order

        Parameters
        ----------
        n : Zernike order
        m : Zernike sub-order
        """
        if m != 0:
            return (2 * (n + 1))**(1 / 2)
        else:
            return (n + 1)**(1 / 2)
    

    def r_nm(self, rho, n, m):
        """Calculates an intermediate function needed for the Zernike polynomials
 
        Parameters
        ----------
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
        
        Parameters
        ----------
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
    """class for the Zernike basis vectors, which are orthogonal over the unit disk"""
    # i dont wanna make this rn imma do it later :D


class KBasis(ZernikeParent):
    """class for the Zernike-like K basis vectors, which are orthogonal over a regular polygon with p sides and a radius, center to corner distance, of R0"""

    def __init__(self, n, m, p, R0):
        """
        Parameters
        ----------
        n : Zernike order
        m : Zernike sub-order
        p : number of sides on the regular polygon
        R0 : radius of the regular polygon
        """
        # inheriting from ZernikeBase
        super().__init__(n, m)

        self.p = p
        self.R0 = R0
        self.alpha = np.pi / p # 2 * alpha is the angle spanned by a single sector of the disk/polyon
        self.side_length = R0 / 2**(1/2) # length of one side of the regular polygon

    def show(self):
        print("n:", self.n, "\nm:", self.m, "\np", self.p, "\nR0:", self.R0)

    # mapping functions
    def u_alpha(self, var1, var2 = np.nan):
        """Generates the intermediate function U_alpha from theta (polygon polar) if a single argument is passed or x,y (polygon cartesian) if two arguments are passed. This function is a piece-wise linear function with slope = 1 except at the each theta where a corner occurs. Immediately before the corner, the value of U_alpha is at a maximum of alpha and, immediately after the corner, the value of U_alpha is at a minimum of -alpha.
    
    Parameters
    ----------
    var1 :
        one argument passed - theta : angular polar coordinate on the regular polygon
        two arguments passed - x : x cartesian coordinate on the regular polygon
    var2 :
        one argument passed - nan : not used, only a holder
        two agrummets passed - y : y cartesian coodinate on the regular polygon
    """
        if var2 == np.nan:
            theta = var1
            return self.u_alpha_theta(theta)
        else:
          x = var1
          y = var2
          return self.u_alpha_xy(x,y)

    def u_alpha_theta(self, theta):
        """Generates the intermediate function U_alpha from theta (polygon polar). This function is a piece-wise linear function with slope = 1 except at the each theta where a corner occurs. Immediately before the corner, the value of U_alpha is at a maximum of alpha and, immediately after the corner, the value of U_alpha is at a minimum of -alpha.

        Parameters
        ----------
        theta : angular polar coordinate on the regular polygon
        """
        drop = (theta + self.alpha) / (2 * self.alpha)
        u_alpha = theta - drop.astype(int) * (2 * self.alpha)
        
        return u_alpha


    def u_alpha_xy(self, x, y):
        """Generates the intermediate function U_alpha from x,y (polygon cartesian). See the description for u_alpha_phi for more info on U_alpha.

        Parameters
        ----------
        x : x cartesian coordinate on the regular polygon
        y : y cartesian coordinate on the regular polygon
        """
        return

    def k_nm(self, x, y):
        """Maps the Zernike polynomial to be orthogonal on a regular polygon.
        Parameters
        ----------
        x : x cartesian coordinate on the regular polygon
        y : y cartesian coordiante on the regular polygon
        """
        # converting from x,y (polygon cartesian) -> r,theta (polygon polar)
        r = (x**2 + y**2)**(1/2)
        theta = np.arctan2(y,x) + (2 * np.pi * (y < 0)) # addition of boolean multiplication is needed to fix mapping done with np.arctan2

        

        print(r,theta, self.n, self.m)



k1 = KBasis(0,0,4,1)
k1.show()
print(k1.u_alpha_theta(np.array([0,0,0])))
