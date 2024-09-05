from . import zernike_like
from .zernike_like import (base_input, ZernikeParent,
                          ZBasis, KBasis,
                          ana_cks, num_cks)
import csv
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as la
import time


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

        return (ck * K).astype(np.float64)


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
            toReturn += self.gen_k_n(i, load_type)

        return toReturn

    def plotter(self, z, name = '', title = ''):
        """Plots the values of z over the regular polygon.

        Parameters
        ----------
        z : np.array
            value of a function at that given x,y
        name : str
            file name of the picture you want to save
        """
        fig, ax = plt.subplots()
        plot = ax.contourf(self.x, self.y, z)
        cbar = fig.colorbar(plot)
        plot.set_clim(-10,6)

        if name != '' :
            plt.title(title)
            plt.savefig("output/images/"+name, dpi = 600)
        #plt.clf()
        plt.close('all')

    def plotter_k_each(self, n, load_type):
        """Plots the value of each K function until the desired
        degree is reached.

        n : int
            Zernike order
        load_type : {"ana", "num"}, str
            chooses the integration scheme the precalculated
            ck_nm will be chosen from
        """
        Fa = base_input(self.x, self.y)
        self.plotter(Fa)
        
        Fk = 0
        for i in range(n+1):
            Fk += self.gen_k_n(i, load_type)
            print("n:", i)
            self.plotter(Fk)

    def err_l2_calc(self, n, load_type):
        """Calculates the relative L2 error of the difference
        between the analytical function and the K
        approximation. Each subsequent error calculation uses
        the combined sum of all previous K_n approximations.
            
        Parameters
        ----------
        n : int
            Zernike Order
        load_type : {"ana", "num"}, str
            chooses the integration scheme the precalculated
            ck_nm will be chosen from
        
        Notes
        -----
        The error calculations are done using relative error
        between the difference of the analytical function and
        the approximation using the various K's as basis
        vectors.

        Inside the calculations, `Fa` is the value of the
        analytical function at the specified x,y, whereas `Fk`
        is the value of the current approximation using K of
        the relevant order as the basis vectors.
        """
        Fa = base_input(self.x, self.y)
        
        l2_errs = []
        Fk = 0

        for i in range(n+1):
            print("n:", i)

            Fk += self.gen_k_n(i, load_type)
            Dk = Fa - Fk

            l2_err = la.norm(Dk, 2) / la.norm(Fa, 2)
            l2_errs.append(l2_err)

        return l2_errs


    def err_linf_calc(self, n, load_type):
        """Calculates the relative Linf error of the difference
        between the analytical function and the K
        approximation. Each subsequent error calculation uses
        the combined sum of all previous K_n approximations.
            
        Parameters
        ----------
        n : int
            Zernike Order
        load_type : {"ana", "num"}, str
            chooses the integration scheme the precalculated
            ck_nm will be chosen from

        Notes
        -----
        The error calculations are done using relative error
        between the difference of the analytical function and
        the approximation using the various K's as basis
        vectors.

        Inside the calculations, `Fa` is the value of the
        analytical function at the specified x,y, whereas `Fk`
        is the value of the current approximation using K of
        the relevant order as the basis vectors.
        """
        Fa = base_input(self.x, self.y)

        linf_errs = []
        Fk = 0 

        for i in range(n+1):
            print("n:", i)

            Fk += self.gen_k_n(i, load_type)
            Dk = Fa - Fk

            linf_err = la.norm(Dk, np.inf) / la.norm(Fa, np.inf)
            linf_errs.append(linf_err)

        return linf_errs


    def kit_and_caboodle(self, n):
        """This function is what I ran to generate all the
        pictures and errors at once. It is an amalgamation of
        most of the above functions, but is created like this,
        so it is easy to run when cloned.
        
        Parameters
        ----------
        n : int
            Zernike order
        
        Notes
        -----
        This function...
            - plots the analytical input
            - plots the K approximation using scipy's ck's
            - plots the K approximation using my ck's
            - gets the l2 error of both K approximations
            - gets the linf error of both K approximation
        """
        
        # initializing
        Fa = base_input(self.x, self.y)
        Fa_l2_norm = la.norm(Fa, 2)
        Fa_linf_norm = la.norm(Fa, np.inf)

        Fk_sci = 0
        Fk_my = 0


        # plotting the analytical input
        print("Plotting the function analytically.")
        self.plotter(Fa, "ana.png", "Analytical")

        print("Starting the Calculations...")
        for i in range(n+1):
            t0 = time.time()
            print("----- n = {} -----".format(i))


            # getting each K for the respective order
            print("calculating approximations")
    
            Fk_sci += self.gen_k_n(i, "ana")
            Fk_my += self.gen_k_n(i, "num")

            tf = time.time()
            print(tf - t0)


            # plotting the approximations
            print("plotting approximations")

            self.plotter(Fk_sci, "k_sci/k_sci{}.png".format(i))
            self.plotter(Fk_my, "k_my/k_my{}.png".format(i))

            tf = time.time()
            print(tf - t0)


            # calculating errors
            print("calculating errors")

            Dk_sci = Fa - Fk_sci
            Dk_my = Fa - Fk_my

            # l2 errors
            l2_err_sci = la.norm(Dk_sci, 2) / Fa_l2_norm
            l2_err_my = la.norm(Dk_my, 2) / Fa_l2_norm

            # linf errors
            linf_err_sci = (la.norm(Dk_sci, np.inf)
                            / Fa_linf_norm)    
            linf_err_my = (la.norm(Dk_my, np.inf)
                           / Fa_linf_norm)
            ''' 
            # writing to output files
            f = open("output/errors/l2_sci.txt", 'a')
            f.write(str(l2_err_sci) + '\n')
            f.close

            f = open("output/errors/l2_my.txt", 'a')
            f.write(str(l2_err_my) + '\n')
            f.close

            f = open("output/errors/linf_sci.txt", 'a')
            f.write(str(linf_err_sci) + '\n')
            f.close

            f = open("output/errors/linf_my.txt", 'a')
            f.write(str(linf_err_my) + '\n')
            f.close
            '''

            tf = time.time()
            print(tf - t0)
