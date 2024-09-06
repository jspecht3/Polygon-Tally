import numpy as np
from math import factorial
from scipy.integrate import nquad
import csv

# base input function
def base_input(x, y):
    """Base input used for the approximation. The basis
    vectors for the subsequent classes will be used to
    decompose this function into a sum of coefficeints
    multiplied by the corresponding basis vector.

    Parameters
    ----------
    x : float
        x cartesian coordinate
    y : float
        y cartesian coordinate

    Notes
    -----
    This function is completely arbitary, I chose it on a
    whim and have no justification as to why I chose it
    other than it is a polynomial with a non-trivial amount
    of complexity.
    """
    return (2 * x**2) - (y**2) + (x**2 * y) - (4 * x * y**2) + (5 * x * y) - (3 * x) + (5 * y)


# precalculated ck values
def parse_ck_csv():
    """Parses data/ck_num_vs_ana.txt and creates a dictionary
    containing those values to shorten runtimes."""
    # lists of each data type
    ck_n = []
    ck_m = []
    num_ck = []
    ana_ck = []

    # filling the lists
    current_path = __file__
    split = current_path.split("zernike_like.py")[0]
    data_path = split + "data/ck_num_vs_ana.txt"

    with open(data_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            ck_n.append(row[0])
            ck_m.append(row[1])
            ana_ck.append(row[2])
            num_ck.append(row[3])

    # making the empty dictionaries
    ana_ck_dict = {}
    num_ck_dict = {}

    for n in range(60):
        ana_ck_m_dict = {}
        num_ck_m_dict = {}

        for m in np.arange(-n,n+1,2):
            ana_ck_m_dict[m] = 0
            num_ck_m_dict[m] = 0

        ana_ck_dict[n] = ana_ck_m_dict
        num_ck_dict[n] = num_ck_m_dict

    # filling the dictionaries with the right values
    for i in range(len(ck_n)):
        n = int(ck_n[i])
        m = int(ck_m[i])

        ana_ck_dict[n][m] = float(ana_ck[i])
        num_ck_dict[n][m] = float(num_ck[i])

    return ana_ck_dict, num_ck_dict

ana_cks, num_cks = parse_ck_csv()


class ZernikeParent:
    """Parent class for the ZBasis, KBasis, HBasis classes""" 

    def __init__(self, n, m):
        """
        Parameters
        ----------
        n : int
            Zernike order
        m : int
            Zernike sub-order
        """
        # checking validity
        self.n_checker(n)
        self.m_checker(n, m)
        
        # initalizing parameters
        self.n = n
        self.m = m


    # value checkers 
    def n_checker(self, n):
        """Verfies the n associated with the Zernike order is
        valid.

        Parameters
        ----------
        n : int
            Zernike order
        """
        if type(n) is not int:
            raise TypeError(
                    "invalid n; make sure n is an integer"
                    )
        
        if n < 0:
            raise ValueError(
                    "invalid n; make sure n >= 0"
                    )
    
        return


    def m_checker(self, n, m):
        """Verifies the m, Zernike sub-order, associated with
        the given n, Zernike order, is valid.
        
        Parameters
        ----------
        n : int
            Zernike order
        m : int
            Zernike sub-order
        """
        if m not in np.arange(-n, n+1, 2):
            raise ValueError(
                    "invalid m; enter an m from [-n, -n+2, -n+4, ... , n]"
                    )

        return


    def rho_checker(self, rho):
        """Verifies the given rho is valid for the domain of the
        Zernike polynomials. rho ranges from [0,1].

        Parameters
        ----------
        rho : float
            radial polar coordinate on the unit disk
        """
        if rho > 1.0 or rho < 0.0:
            raise ValueError(
                    "rho out of bounds; rho ranges from [0,1]"
                    )

        return
        
    def phi_checker(self, phi):
        """Verifies the given phi is valid for the domain of the
        Zernike polynomials. phi ranges from [0, 2*pi).

        Parameters
        ----------
        phi : float
            angular polar coordinate on the unit disk
        """
        if phi > 2 * np.pi or phi < 0.0:
            raise ValueError(
                    "phi out of bounds; phi ranges from [0,2pi)"
                    )
    
        return
    

    # zernike intermediate functions
    def n_nm(self, n, m):
        """Calculates the normalization constant for the Zernike
        polynomials.

        Parameters
        ----------
        n : int
            Zernike order
        m : int
            Zernike sub-order
        """
        if m != 0:
            return (2 * (n + 1))**(1 / 2)
        else:
            return (n + 1)**(1 / 2)
    

    def r_nm(self, rho, n, m):
        """Calculates an intermediate function needed for the
        Zernike polynomials.
 
        Parameters
        ----------
        rho : float
            radial polar cordinate on the unit disk
        n : int
            Zernike order
        m : int
            Zernike sub-order
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
        """Creates the zernike function associated with the
        respective n, Zernike order, and m, Zernike sub-order,
        values
        
        Parameters
        ----------
        rho : float
            radial polar coordinate on the unit disk
        phi : float
            angular polar coordinate on the unit disk
        n : int
            Zernike order
        m : int
            Zernike sub-order
        """
        # checking
        #self.rho_checker(rho)
        #self.phi_checker(phi)

        # value of the Zernike polynomial
        if m >= 0:
            return self.n_nm(n, m) * self.r_nm(rho, n, m) * np.cos(m * phi)
        else:
            return -1 * self.n_nm(n, m) * self.r_nm(rho, n, m) * np.sin(m * phi)
    

class ZBasis(ZernikeParent):
    """Class for the Zernike basis vectors, which are orthogonal
    over the unit disk"""
    # i dont wanna make this rn imma do it later


class KBasis(ZernikeParent):
    """Class for the Zernike-like K basis vectors, which are
    orthogonal over a regular polygon with p sides and a
    radius, center to corner distance, of R0"""

    def __init__(self, n, m, num_sides, polygon_radius):
        """
        Parameters
        ----------
        n : int
            Zernike order
        m : int
            Zernike sub-order
        num_sides : int
            number of sides on the regular polygon
        polygon_radius : float
            radius of the regular polygon
        
        Notes
        -----
        alpha : float
            half the angle spanned by a single sector of the
            disk/polygon. As there are the same number of sides
            as there are sectors and alpha spans half a sector,
            2*alpha can be found by dividing 2*pi by num_sides.
        side_length : float
            side length of a side of the regular polygon 
        ck : float
            the factorization contant used in the decomposition
            of a function into a series of KBasis vectors
        ck_type : str
            ck can be calculated in a few ways, so this shows
            the current way it was calculated
        """
        # inheriting from ZernikeBase
        super().__init__(n, m)

        self.num_sides = num_sides
        self.polygon_radius = polygon_radius
        self.alpha = np.pi / num_sides
        self.side_length = polygon_radius / 2**(1/2)
        self.ck = None
        self.ck_type = None
        self.mesh_size = None

    def show(self):
        print("n:", self.n)
        print("m:", self.m)
        print("Number of Sides:", self.num_sides)
        print("Polygon Radius:", self.polygon_radius)
        print("ck value:", self.ck)
        print("ck Calculation Method:", self.ck_type)
        print("Mesh Size:", self.mesh_size)


    # start of coordinate transforms
    def cart_to_polar(self, x, y):
        """Takes the cartesian coordinates on the polygon (x,y)
        and returns the polar coordinates on the polygon
        (r,theta).
        
        Parameters
        ----------
        x : float
            x cartesian coordinate on the regular polygon
        y : float
            y cartesian coordinate on the regular polygon

        Notes
        -----
        The addition of the boolean multiplication when
        calculating theta is needed to fix the mapping done
        with np.arctan2, which is unable to determine the
        quadrant of the coordinates correctly.
        """
        r = (x**2 + y**2)**(1/2)
        theta = np.arctan2(y,x) + (2 * np.pi * (y < 0))
        
        return r, theta


    def polar_to_cart(self, r, theta):
        """Takes the polar coordinates on the polygon (r,theta)
        and returns the cartesian coorinates on the polygon
        (x,y)

        Parameters
        ----------
        r : float
            radial polar coordinates on the regular polygon
        theta : float
            angular polar coordinates on the regular polygon
        """
        x = r * np.cos(theta)
        y = r * np.sin(theta)

        return x, y
    # end of coordinate transforms


    # start of variable radius functions
    def u_alpha_theta(self, theta):
        """Generates the intermediate function U_alpha from
        theta (polygon polar). This function is a piece-wise
        linear function with slope = 1 except at the each theta 
        where a corner occurs. Immediately before the corner,
        the value of U_alpha is at a maximum of alpha and,
        immediately after the corner, the value of U_alpha is
        at a minimum of -alpha.

        Parameters
        ----------
        theta : float
            angular polar coordinate on the regular polygon
        
        Notes
        -----
        The drop variable is used to find the locations of the
        corners, which are also the points of discontinuity, of
        the polygon.
        """
        drop = (theta + self.alpha) / (2 * self.alpha)
        u_alpha = theta - drop.astype(int) * (2 * self.alpha) 
#        u_alpha = theta - int(drop) * (2 * self.alpha)
        
        return u_alpha
    

    def u_alpha_xy(self, x, y):
        """Generates the intermediate function U_alpha from x,y 
        (polygon cartesian).

        Parameters
        ----------
        x : float
            x cartesian coordinate on the regular polygon
        y : float
            y cartesian coordinate on the regular polygon
        
        See Also
        --------
        u_alpha_theta : same function in polar coordinates
        """
        _, theta = self.cart_to_polar(x, y)
        u_alpha = self.u_alpha_theta(theta)

        return u_alpha

    
    def variable_radius_theta(self, theta):
        """Finds the variable radius on the regular polygon for
        a given theta. The variable radius is the distance from
        the center of the polygon to the edge of the polygon at
        that given theta.

        Parameters
        ----------
        theta : float
            angular polar coordinate on the regular polygon
        """
        numerator = self.polygon_radius * np.cos(self.alpha)
        denominator = np.cos(self.u_alpha_theta(theta))
        
        variable_radius = numerator / denominator

        return variable_radius


    def variable_radius_xy(self, x, y):
        """Finds the variable radius on the regular polygon for
        a given x,y.

        Parameters
        ----------
        x : float
            x cartesian coordinate on the regular polygon
        y : float
            y cartesian coordinate on the regular polygon

        See Also
        --------
        variable_radius_theta : same function for polar
        coordinates
        """
        _, theta = self.cart_to_polar(x, y)
        variable_radius = self.variable_radius_theta(theta)

        return variable_radius
    # end of variable radius functions


    # start of coordinate generators
    ## start coordinate generators
    def _check_mesh_size_exists(self):
        """Verifies the KBasis has a mesh size associated for
        numerical integration.
        """
        if self.mesh_size == None:
            raise AttributeError("the attribute 'self.mesh_size' is undefined.")

    
    def _gen_rhophi(self):
        """Generates a matrix of mesh_size by mesh_size over
        the polar coordinate space of the disk.
        """
        self._check_mesh_size_exists()

        rho = np.linspace(0, 1, self.mesh_size)
        phi = np.linspace(0, 2 * np.pi, self.mesh_size)
        
        rho, phi = np.meshgrid(rho, phi)
        
        self._rho = rho
        self._phi = phi


    def _gen_uv(self):
        """Generates a matrix of mesh_size by mesh_size over
        the cartesian coordinate space of the disk.
        """
        self._check_mesh_size_exists()

        if not hasattr(self, "_rho"):
            raise AttributeError("No rho, phi defined")
        
        self._u = self._rho * np.cos(self._phi)
        self._v = self._rho * np.sin(self._phi)

    
    def _gen_rtheta(self):
        """Generates a matrix of mesh_size by mesh_size over
        the polar coordinate space of the polygon.
        """
        self._check_mesh_size_exists()

        if not hasattr(self, "_rho"):
            raise AttributeError("No rho, phi defined")

        var_rad = self.variable_radius_theta(self._phi)
        self._r = var_rad * self._rho
        self._theta = self._phi


    def _gen_xy(self):
        """Generates a matrix of mesh_size by mesh_size over
        the cartesian coordinate space of the polygon.
        """ 
        self._check_mesh_size_exists()

        if not hasattr(self, "_r"):
            raise AttributeError("No r, theta defined")

        self._x = self._r * np.cos(self._theta)
        self._y = self._r * np.sin(self._theta)
    ## end coordinate generators


    ## start polar differentials
    def _gen_polar_diffs_disk(self):
        """Generates the polar differentials for the disk"""
        if not hasattr(self, "_rho"):
            raise AttributeError("No rho, phi defined")

        drho = self._rho[:,1]
        _, self._drho = np.meshgrid(drho, drho)
        self._dphi = self._phi[1,1]
        
    
    def _gen_polar_diffs_polygon(self):
        """Generates the polar differentials for the polygon"""

        if not hasattr(self, "_r"):
            raise AttributeError("No r, theta defined")

        dr = self._r[:,1]
        _, self._dr = np.meshgrid(dr, dr)
        self._dtheta = self._theta[1,1]
    ## end polar differentials


    def _gen_all_coords(self, mesh_size):
        """Generates polar, cartesian coordinates for the
        disk, polygon.
        """
        self.mesh_size = mesh_size
        self._gen_rhophi()
        self._gen_uv()
        self._gen_rtheta()
        self._gen_xy()
   

    def _gen_all_diffs(self, mesh_size):
        """Generates all polar differentials used in the
        numerical integration for the disk, polygon."""
        self.mesh_size = mesh_size
        self._gen_polar_diffs_disk()
        self._gen_polar_diffs_polygon()

    def _gen_all_cd(self, mesh_size):
        "Generates all coordinates and differentials"""
        self._gen_all_coords(mesh_size)
        self._gen_all_diffs(mesh_size)


    def check_coords(self):
        """Checks which coordinates are defined as attributes.
        """
        print("------ Coordinates ------")
        if hasattr(self, "_rho"):
            print("Polar Disk :      DEFINED")
        else:
            print("Polar Disk :    UNDEFINED")

        if hasattr(self, "_u"):
            print("Cart Disk  :      DEFINED")
        else:
            print("Polar Disk :    UNDEFINED")

        if hasattr(self, "_r"):
            print("Polar Poly :      DEFINED")
        else:
            print("Polar Poly :    UNDEFINED")

        if hasattr(self, "_x"):
            print("Cart Poly  :      DEFINED")
        else:
            print("Cart Poly  :    UNDEFINED")
        
        print("-- Polar Differentials --") 
        if hasattr(self, "_drho"):
            print("Diffs Disk :      DEFINED")
        else:
            print("Diffs Disk :    UNDEFINED")
 
        if hasattr(self, "_dr"):
            print("Diffs Poly :      DEFINED")
        else:
            print("Diffs Poly :    UNDEFINED")
    # end of coordinate generators


    # start of basis vectors
    def k_nm(self, x, y):
        """Maps the Zernike polynomial to be orthogonal on a
        regular polygon.
        
        Parameters
        ----------
        x : float
            x cartesian coordinate on the regular polygon
        y : float
            y cartesian coordiante on the regular polygon
        
        Notes
        -----
        rho : float
            radial polar coordinate on the unit disk
            It should be noted that k_nm is simply zernike_nm,
            but mapped to the polygon. K_nm is internally
            mapping the coordinates of the polygon to the
            corresponding coordinates on the disk.
        """
        r, theta = self.cart_to_polar(x, y)
        variable_radius = self.variable_radius_theta(theta)

        rho = r / variable_radius
        phi = theta
        
        k_nm_value = self.zernike_nm(rho, phi, self.n, self.m)

        return k_nm_value

    
    def num_ck_nm(self, mesh_size, function = base_input):
        """Numerically calculates the coefficient, ck, for the
        respective k_nm using a simple, cenetered, Riemann
        numerical integration scheme.

        Parameters
        ---------- 
        mesh_size : int
            size of the mesh used for the numerical integration
        function : func(x,y)
            function to be approximated by decomposition into a
            sum of ck's times the k_nm's.
            function = Sigma_(n=0)^(n=infty) [ck_nm * k_nm]
        """
        dmu = self._rho * self._drho * self._dphi / np.pi

        ck = np.sum(dmu * function(self._x, self._y) * self.k_nm(self._x, self._y))

        self.ck = ck
        self.ck_type = "Numerical Integration w/ Riemann Sum"

        return ck


    def ana_ck_nm(self, function = base_input):
        """Analytically calculates the coefficient, ck, for the
        respective k_nm using scipy.

        Parameters
        ----------
        function : func(x,y)
            function to be approximated by decomposition into a
            sum of ck's times the k_nm's.
            function = Sigma_(n=0)^(n=infty) [ck_nm * k_nm]

        See Also
        --------
        add this once you come up with a name for the function

        Notes
        -----
        THE DESCRIPTION IS WRONG!!! THIS IS NOT AN ANALYTICAL
        INTEGRATION SCHEME. I WANT TO DO OTHER THINGS ATM, SO
        I AM NOT CHANGING THE NAMING CONVENTIONS, BUT THIS IS
        USING SCIPY'S QUADRATURE PACK, WHICH COMES FROM THE
        FORTRAN QUADRATURE RULES.

        THIS ALSO ONLY WORKS WITH A SQUARE FOR NOW. MAYBE
        FOREVER BECAUSE THIS IS ONLY MEANT TO BE A COMPARISON
        AGAINST A RIEMANN SUM.

        The calculation for this can get rather computationally
        expensive, which is especially true when considering
        cumulative, all m's for a given n, approximations at
        high orders. For a given Zernike order n, the there are
        (n+1) associated m values, meaning the runtime per
        cumulative Zernike order increases not only due to the
        increased complexity of higher order Zernikes, but also
        due to the fact there are more polynomials to work with
        in each subsequent order.
        """
        def _toInt(x, y):
            K = self.k_nm(x,y)
            F = function(x,y)

            R = self.variable_radius_xy(x,y)

            dmu = 1 / np.pi / R**2
            
            return K * F * dmu

        sl = self.side_length
        ck = nquad(_toInt, [[-sl,sl],[-sl,sl]])[0]
    
        self.ck = ck
        self.ck_type = "Calculated with scipy"

        return ck


    def load_ck_nm(self, load_type):
        """Loads precalculated ck_values to reduce the runtime
        of the approximations. These values were calculated
        with the same functions that are used in this class.
        (THEN ACTAULLY DO THAT AND MAKE THE NUMERICAL ONE TOO)

        Parameters
        ----------
        load_type : {"ana", "num"}, str
            chooses which integration scheme the precalculated
            ck_nm will be chosen from
        """
        if load_type == "ana":
            self.ck = ana_cks[self.n][self.m]
            self.ck_type = "Analytically Fetched"
        
        elif load_type == "num":
            self.ck = num_cks[self.n][self.m]
            self.ck_type = "Numerically Fetched"

        else:
            raise ValueError(
            "load_type can be either 'ana' or 'num', which loads the analyical or numerical ck, respectively."
            )
        
        return 
"""
    # wrong value, but compiles
    def testing1(self):
        def _toInt(rho, phi):
            theta = phi
            r = rho * self.variable_radius_theta(phi)

            x = r * np.cos(theta)
            y = r * np.sin(theta)

            K = self.k_nm(x, y)
            F = base_input(x, y)
            R = self.variable_radius_xy(x, y)

            dmu = rho / np.pi / R**2

            return K * F * dmu

        ck = nquad(_toInt, [[0 , 1], [0, 2 * np.pi]])[0]
        return ck

    # works
    def testing2(self):
        def _toInt(x, y):
            K = self.k_nm(x, y)
            F = base_input(x, y)

# output for testing
k1 = KBasis            R = self.variable_radius_xy(x, y)

            dmu = 1 / np.pi / R**2
            return K * F * dmu
        
        sl = self.side_length
        ck = nquad(_toInt, [[-sl, sl], [-sl, sl]])[0]

        return ck

    # now compiles
    def testing3(self):
        def _toInt(u, v):
            rho, phi = self.cart_to_polar(u, v) 
            R = self.variable_radius_theta(phi)

            x = u * R
            y = v * R

            K = self.k_nm(x, y)
            F = base_input(x, y)
            dmu = 1 / np.pi / R**2

            return K * dmu * F

        def ub():
            return [-1,1]
        def vb(u):
            a = (1 - u**2)**(1/2)
            return [-a, a]

        ck = nquad(_toInt, [vb,ub])[0]
        return ck
"""


"""
# output for testing
k1 = KBasis(0,0,4,1)
k1.show()
print(k1.testing1())
print(k1.testing2())
print(k1.testing3())
"""
