import openmc
import numpy as np
import matplotlib.pyplot as plt
from hex_rings import radius

# dimesnsions
pin_pitch = 1.180 * 5.5685 / 10
polygon_radius = pin_pitch * 10 * 3**(0.5)

# statepoint
sp = openmc.StatePoint('statepoint.100000.h5')

zern_tally = sp.tallies[1]
poly_tally = sp.tallies[2]


def get_stop(order):
    stops = []
    count = 0
    for i in range(16):
        count += i + 1
        stops.append(count)
    return stops[order]


# zernike fission info
zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_fission_df = zern_fission.get_pandas_dataframe()
zern_fission_n = zern_fission_df['mean']


def zern_fission_me(order):
    stop = get_stop(order)
    ns = zern_fission_n[:stop]
    return openmc.Zernike(ns, radius=1)


# zernike flux info
zern_flux = zern_tally.get_slice(scores=['flux'])
zern_flux_df = zern_flux.get_pandas_dataframe()
zern_flux_n = zern_flux_df['mean']
zern_flux_zs = openmc.Zernike(zern_flux_n, radius=1)


def zern_flux_me(order):
    stop = get_stop(order)
    ns = zern_flux_n[:stop]
    return openmc.Zernike(ns, radius=1)


# polygon fission info
poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
poly_fission_df = poly_fission.get_pandas_dataframe()
poly_fission_n = poly_fission_df['mean']
poly_fission_zs = openmc.Zernike(poly_fission_n, radius=1)


def poly_fission_me(order):
    stop = get_stop(order)
    ns = poly_fission_n[:stop]
    return openmc.Zernike(ns, radius=1)


# polygon flux info
poly_flux = poly_tally.get_slice(scores=['flux'])
poly_flux_df = poly_flux.get_pandas_dataframe()
poly_flux_n = poly_flux_df['mean']
poly_flux_zs = openmc.Zernike(poly_flux_n, radius=1)


def poly_flux_me(order):
    stop = get_stop(order)
    ns = poly_flux_n[:stop]
    return openmc.Zernike(ns, radius=1)


# functions for zern tally
def zern_fission_xy(x, y, order):
    r = (x**2 + y**2)**(1/2) / polygon_radius
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))
    return np.array([zern_fission_me(order)(r, theta)], dtype=float) * np.pi


def zern_flux_xy(x, y, order):
    r = (x**2 + y**2)**(1/2) / polygon_radius
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))
    return np.array([zern_flux_me(order)(r, theta)], dtype=float) * np.pi


# functions for poly tally
num_sides = 6
alpha = np.pi / num_sides


def r_alpha(theta):
    drop = (theta + alpha) / (2 * alpha)
    u_alpha = theta - drop.astype(int) * (2 * alpha)
    return polygon_radius * np.cos(alpha) / np.cos(u_alpha)


def poly_fission_xy(x, y, order):
    r = (x**2 + y**2)**(1/2)
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))

    var_radius = r_alpha(theta)
    rho = r / var_radius
    return np.array([poly_fission_me(order)(rho, theta)], dtype=float) * np.pi


def poly_flux_xy(x, y, order):
    r = (x**2 + y**2)**(1/2)
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))

    var_radius = r_alpha(theta)
    rho = r / var_radius
    return np.array([poly_flux_me(order)(rho, theta)], dtype=float) * np.pi
