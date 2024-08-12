from zernike_like import (base_input, ZernikeParent,
                          ZBasis, KBasis)
import matplotlib.pyplot as plt
import numpy as np


class ApproxParent():
    """Parent Class for the ZApprox, KApprox, and HApprox
    classes"""

    def __init__(self):
        return


class KApprox(ApproxParent):

    def __init__(self, num_sides, polygon_radius, mesh_size):
        """
        Parameters
        ----------
        num_sides : int
            number of sides on the regular polygon
        polygon_radius : float
            radius of the polygon
        mesh_side : int
            the number of elements used in the array used to
            generate the value of each basis vector
        """
        # inheriting from ApproxParent
        super().__init__()

        # simple initializers
        self.num_sides = num_sides
        self.polygon_radius = polygon_radius
        self.mesh_size = mesh_size

        # coordinate initializers
        self.k00 = KBasis(0, 0, self.num_sides, self.polygon_radius)
        self.gen_coordinates(self.mesh_size)


    def gen_coordinates(self, mesh_size):
        """Generates the various coordinates needed for
        subsequent calculations.
        """

        def _gen_rhophi(mesh_size):
            """Generates the polar coordinates on the unit disk
            """

            rho = np.linspace(0, 1, mesh_size)
            phi = np.linspace(0, 2 * np.pi, mesh_size)
            
            rho, phi = np.meshgrid(rho, phi)

            self.rho = rho
            self.phi = phi
            
            return

        def _gen_rtheta(mesh_size):
            """Generates the polar coordinates on the regular
            polygon.
            """
            theta = self.phi

            variable_radius = self.k00.variable_radius_theta(theta)
            r = variable_radius * self.rho

            self.theta = theta
            self.r = r

            return

        def _gen_xy(mesh_size):
            """Generate the cartesian coordinates on the
            regular polygon.
            """
            x = self.r * np.cos(self.theta)
            y = self.r * np.sin(self.theta)

            self.x = x
            self.y = y

            return
 
        _gen_rhophi(mesh_size)
        _gen_rtheta(mesh_size)
        _gen_xy(mesh_size)

        return


    def gen_k_nm(self, n, m, load_type):
        """Generates the values of a K basis vector given a
        specific n and m value.

        Parameters
        ----------
        n : int
            Zernike order
        m : int
            Zernike sub-order
        load_type : {"ana", "num"}, str
            chooses the integration scheme the precalculated
            ck_nm will be chosen from
        """
        basis = KBasis(n, m, self.num_sides, self.polygon_radius)
        basis.load_ck_nm(load_type)

        K = basis.k_nm(self.x, self.y)
        ck = basis.ck

        return ck * K


    def gen_k_n(self, n, load_type):
        """Generates the summed values of all K basis vectors
        given a specific n value.

        Parameters
        ----------
        n : int
            Zernike order
        load_type : {"ana", "num"}, str
            chooses the integration scheme the precalculated
            ck_nm will be chosen from
        """

        toReturn = 0

        for m in np.arange(-n, n+1, 2):
            toReturn += self.gen_k_nm(n, m, load_type)

        return toReturn


    def gen_k_total(self, n, load_type):
        """Generates the approximation of the inputted function
        using the transformed Zernike polynomials, K, as the
        basis vectors

        Parameters
        ----------
        n : int
            Zernike order
        load_type : {"ana", "num"}, str
            chooses the integration scheme the precalculated
            ck_nm will be chosen from
        """
        toReturn = 0

        for i in range(n + 1):
            toReturn += self.gen_k_n(n, load_type)

        return toReturn

    def plotter(self, z):
        fig, ax = plt.subplots()
        plot = ax.contourf(self.x, self.y, z)
        cbar = fig.colorbar(plot)
        plt.savefig("testing.png")
        plt.show()






a = KApprox(4, 100, 1000)
z1 = a.gen_k_n(5,"ana")
