import openmc
import numpy as np
import matplotlib.pyplot as plt
from hex_rings import radius

# dimesnsions
pin_pitch = 1.180 * 5.5685 / 10
polygon_radius = pin_pitch * 10 * 3**(0.5)

# statepoint
sp = openmc.StatePoint('statepoint.10000.h5')

zern_tally = sp.tallies[2]
poly_tally = sp.tallies[3]

# zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
# zern_flux = zern_tally.get_slice(scores=['flux'])

# poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
# poly_flux = poly_tally.get_slice(scores=['flux'])


# getting zernike fission info
zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_fission_df = zern_fission.get_pandas_dataframe()
zern_fission_n = zern_fission_df['mean']
zern_fission_zs = openmc.Zernike(zern_fission_n, radius=1)

# getting polygon fission info
poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
poly_fission_df = poly_fission.get_pandas_dataframe()
poly_fission_n = poly_fission_df['mean']
poly_fission_zs = openmc.Zernike(poly_fission_n, radius=1)


# functions for zern tally
def zern_fission_xy(x, y):
    r = (x**2 + y**2)**(1/2) / polygon_radius
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))
    return zern_fission_zs(r, theta)


# functions for poly tally
num_sides = 6
alpha = np.pi / num_sides


def r_alpha(theta):
    drop = (theta + alpha) / (2 * alpha)
    u_alpha = theta - drop.astype(int) * (2 * alpha)
    return polygon_radius * np.cos(alpha) / np.cos(u_alpha)


def poly_fission_xy(x, y):
    r = (x**2 + y**2)**(1/2)
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))

    var_radius = r_alpha(theta)
    rho = r / var_radius
    return poly_fission_zs(rho, theta)
