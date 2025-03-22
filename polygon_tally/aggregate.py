from . import zernike_like
from .zernike_like import (base_input, ZernikeParent,
                          ZBasis, KBasis,
                          ana_cks, num_cks)
import csv
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
import numpy as np
import numpy.linalg as la
import time
import os

class ApproxParent():
    """Parent Class for the ZApprox, KApprox, and HApprox
    classes"""

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
        # simple initializers
        self.num_sides = num_sides
        self.polygon_radius = polygon_radius
        self.mesh_size = mesh_size
        return

    def plotter(self, z, bounds = [], title = '', name = ''):
        """Plots the values of z over the regular polygon.

        Parameters
        ----------
        z : np.array
            value of a function at that given x,y
        title : str
            figure title
        name : str
            file name of the picture you want to save
        """
        cmap = mpl.colormaps['viridis']
        norm = mpl.colors.Normalize(bounds[0], bounds[1])
        sm = ScalarMappable(norm=norm, cmap=cmap)
        sm.set_array([])

        fig, ax = plt.subplots()
        plot = ax.contourf(self.x, self.y, z, cmap=cmap, levels=1000)
        cbar = fig.colorbar(sm, ax=ax)
        
#        if bounds != []:
#            plot.set_clim(bounds)
       
        graph_lim = 1.1 * self.polygon_radius
        ax.set_xlim(-graph_lim, graph_lim)
        ax.set_ylim(-graph_lim, graph_lim)
        ax.set_title(title)
        
        if name != '' :
            plt.savefig(name, dpi = 600)
        #plt.close('all')
        plt.close()

class ZApprox(ApproxParent):

    def __init__(self, num_sides, polygon_radius, mesh_size):
        """
        Parameters
        ----------
        num_sides : int
            number of sides on the regular polygon
        polygon_radius : float
            radius of the polygon
        mesh_size : int
            the number of elements used in the array used to
            generate the value of each basis vector
        """
        # inheriting from ApproxParent
        super().__init__(num_sides, polygon_radius, mesh_size)

        # coordinate initializers
        self.z00 = ZBasis(0, 0, self.num_sides, self.polygon_radius)
        self.z00._gen_all_coords(mesh_size)
        self.rho, self.phi = self.z00._rho, self.z00._phi
        self.x, self.y = self.z00._x, self.z00._y
        return

    def gen_z_nm(self, n, m, cz_method, cz_type, function = base_input):
        """Generates the values of a Z basis vector given a
        specific n and m value.

        Parameters
        ----------
        n : int
            Zernike order
        m : int
            Zernike sub-order
        cz_method : {"calc", "load"}, str
            chooses to either calculate or load the zk
        cz_type : {"num", "ana"}, str
            chooses the calculation scheme for the zk
        function : func(x,y)
            the function to be expanded with the ZBasis
        """
        basis = ZBasis(n, m, self.num_sides, self.polygon_radius)
        basis._gen_all_cd(self.mesh_size)

        # add support for other methods here
        if (cz_method == "calc") and (cz_type == "num"):
            cz = basis.num_cz_nm(self.mesh_size)
        else:
            raise ValueError("Support isn't here for that rn")
        
        return (cz * basis.z_nm(self.x, self.y))
    
    def gen_z_n(self, n, cz_method, cz_type, function = base_input):
        """Generates the summed values of all Z basis vectors
        given a specific n value.

        Parameters
        ----------
        n : int
            Zernike order
        cz_method : {"calc", "load"}, str
            chooses to either calculate or load the cz
        cz_type : {"num", "ana"}, str
            chooses the calculation scheme for the cz
        function : func(x,y)
            the function to be expanded with the ZBasis
        """
        toReturn = 0

        for m in np.arange(-n, n+1, 2):
            toReturn += self.gen_z_nm(n, m, cz_method, cz_type, function)
        return toReturn

    def gen_z_total(self, n, cz_method, cz_type, function = base_input):
        """Generates the approximation of the inputted function
        using the Zernike polynomials as the basis vectors

        Parameters
        ----------
        n : int
            Zernike order
        cz_method : {"calc", "load"}, str
            chooses to either calculate or load the cz
        cz_type : {"num", "ana"}, str
            chooses the calculation scheme for the cz
        function : func(x,y)
            the function to be expanded with the ZBasis
        """
        toReturn = 0

        for i in range(n + 1):
            toReturn += self.gen_z_n(i, cz_method, cz_type, function)
        return toReturn

    def err_l2_calc(self, n, cz_method, cz_type, function = base_input):
        """Calculates the relative L2 error of the difference
        between the analytical function and the Z
        approximation. Each subsequent error calculation uses
        the combined sum of all previous Z_n approximations.
            
        Parameters
        ----------
        n : int
            Zernike Order
        cz_method : {"calc", "load"}, str
            chooses to either calculate or load the cz
        cz_type : {"num", "ana"}, str
            chooses the calculation scheme for the cz
        function : func(x,y)
            the function to be expanded with the ZBasis
        
        Notes
        -----
        The error calculations are done using relative error
        between the difference of the analytical function and
        the approximation using the various Z's as basis
        vectors.

        Inside the calculations, `Fa` is the value of the
        analytical function at the specified x,y, whereas `Fz`
        is the value of the current approximation using K of
        the relevant order as the basis vectors.
        """
        Fa = base_input(self.x, self.y)
        
        l2_errs = []
        Fz = 0

        for i in range(n+1):
            print("n:", i)

            Fz += self.gen_z_n(i, cz_method, cz_type, function)
            Dz = Fa - Fz

            l2_err = la.norm(Dz, 2) / la.norm(Fa, 2)
            l2_errs.append(l2_err)

        return l2_errs

    def err_linf_calc(self, n, cz_method, cz_type, function = base_input):
        """Calculates the relative Linf error of the difference
        between the analytical function and the Z
        approximation. Each subsequent error calculation uses
        the combined sum of all previous Z_n approximations.
            
        Parameters
        ----------
        n : int
            Zernike Order
        cz_method : {"calc", "load"}, str
            chooses to either calculate or load the zk
        cz_type : {"num", "ana"}, str
            chooses the calculation scheme for the zk
        function : func(x,y)
            the function to be expanded with the ZBasis

        Notes
        -----
        The error calculations are done using relative error
        between the difference of the analytical function and
        the approximation using the various Z's as basis
        vectors.

        Inside the calculations, `Fa` is the value of the
        analytical function at the specified x,y, whereas `Fz`
        is the value of the current approximation using Z of
        the relevant order as the basis vectors.
        """
        Fa = base_input(self.x, self.y)

        linf_errs = []
        Fz = 0

        for i in range(n+1):
            print("n:", i)

            Fz += self.gen_z_n(i, cz_method, cz_type, function)
            Dz = Fa - Fz

            linf_err = la.norm(Dz, np.inf) / la.norm(Fa, np.inf)
            linf_errs.append(linf_err)

        return linf_errs

class KApprox(ApproxParent):

    def __init__(self, num_sides, polygon_radius, mesh_size):
        """
        Parameters
        ----------
        num_sides : int
            number of sides on the regular polygon
        polygon_radius : float
            radius of the polygon
        mesh_size : int
            the number of elements used in the array used to
            generate the value of each basis vector
        """
        # inheriting from ApproxParent
        super().__init__(num_sides, polygon_radius, mesh_size)

        # coordinate initializers
        self.k00 = KBasis(0, 0, self.num_sides, self.polygon_radius)
        self.k00._gen_all_coords(mesh_size)
        self.rho, self.phi = self.k00._rho, self.k00._phi
        self.x, self.y = self.k00._x, self.k00._y
        return

    def gen_k_nm(self, n, m, ck_method, ck_type, function = base_input):
        """Generates the values of a Z basis vector given a
        specific n and m value.

        Parameters
        ----------
        n : int
            Zernike order
        m : int
            Zernike sub-order
        ck_method : {"calc", "load"}, str
            chooses to either calculate or load the ck
        ck_type : {"num", "ana"}, str
            chooses the calculation scheme for the ck
        function : func(x,y)
            the function to be expanded with the ZBasis
        """
        basis = KBasis(n, m, self.num_sides, self.polygon_radius)
        basis._gen_all_cd(self.mesh_size)

        # add support for other methods here
        if (ck_method == "calc") and (ck_type == "num"):
            ck = basis.num_ck_nm(self.mesh_size, function)
        else:
            raise ValueError("Support isn't here for that rn")
        
        return (ck * basis.k_nm(self.x, self.y))

    def gen_k_n(self, n, ck_method, ck_type, function = base_input):
        """Generates the summed values of all K basis vectors
        given a specific n value.

        Parameters
        ----------
        n : int
            Zernike order
        ck_method : {"calc", "load"}, str
            chooses to either calculate or load the ck
        ck_type : {"num", "ana"}, str
            chooses the calculation scheme for the ck
        function : func(x,y)
            the function to be expanded with the ZBasis
        """
        toReturn = 0

        for m in np.arange(-n, n+1, 2):
            toReturn += self.gen_k_nm(n, m, ck_method, ck_type, function)
        return toReturn

    def gen_k_total(self, n, ck_method, ck_type, function = base_input):
        """Generates the approximation of the inputted function
        using the transformed Zernike polynomials, K, as the
        basis vectors

        Parameters
        ----------
        n : int
            Zernike order
        ck_method : {"calc", "load"}, str
            chooses to either calculate or load the ck
        ck_type : {"num", "ana"}, str
            chooses the calculation scheme for the ck
        function : func(x,y)
            the function to be expanded with the ZBasis
        """
        toReturn = 0

        for i in range(n + 1):
            toReturn += self.gen_k_n(i, ck_method, ck_type, function)
        return toReturn

    def err_l2_calc(self, n, ck_method, ck_type, function = base_input):
        """Calculates the relative L2 error of the difference
        between the analytical function and the K
        approximation. Each subsequent error calculation uses
        the combined sum of all previous K_n approximations.
            
        Parameters
        ----------
        n : int
            Zernike Order
        ck_method : {"calc", "load"}, str
            chooses to either calculate or load the ck
        ck_type : {"num", "ana"}, str
            chooses the calculation scheme for the ck
        function : func(x,y)
            the function to be expanded with the ZBasis
        
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

            Fk += self.gen_k_n(i, ck_method, ck_type, function)
            Dk = Fa - Fk

            l2_err = la.norm(Dk, 2) / la.norm(Fa, 2)
            l2_errs.append(l2_err)

        return l2_errs

    def err_linf_calc(self, n, ck_method, ck_type, function = base_input):
        """Calculates the relative Linf error of the difference
        between the analytical function and the K
        approximation. Each subsequent error calculation uses
        the combined sum of all previous K_n approximations.
            
        Parameters
        ----------
        n : int
            Zernike Order
        ck_method : {"calc", "load"}, str
            chooses to either calculate or load the ck
        ck_type : {"num", "ana"}, str
            chooses the calculation scheme for the ck
        function : func(x,y)
            the function to be expanded with the ZBasis

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

            Fk += self.gen_k_n(i, ck_method, ck_type, function)
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
