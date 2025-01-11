import openmc
import numpy as np
import matplotlib.pyplot as plt
from hex_rings import radius

import matplotlib as mpl
from matplotlib.cm import ScalarMappable

# dimesnsions
pin_pitch = 1.180 * 5.5685 / 10
polygon_radius = pin_pitch * 10 * 3**(0.5)

# statepoint
sp = openmc.StatePoint('statepoint.100000.h5')

zern_tally = sp.tallies[1]
poly_tally = sp.tallies[2]

zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_flux = zern_tally.get_slice(scores=['flux'])

poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
poly_flux = poly_tally.get_slice(scores=['flux'])


# plotting zernike
zern_fission_df = zern_fission.get_pandas_dataframe()
zern_flux_df = zern_flux.get_pandas_dataframe()

zern_fission_n = zern_fission_df['mean']
zern_flux_n = zern_flux_df['mean']

zern_fission_zs = openmc.Zernike(zern_fission_n, radius=polygon_radius)
zern_flux_zs = openmc.Zernike(zern_flux_n, radius=polygon_radius)

azimuths = np.linspace(0, 2*np.pi, 100)
zeniths = np.linspace(0, polygon_radius, 100)
r, theta = np.meshgrid(zeniths, azimuths)

z_fission = zern_fission_zs(zeniths, azimuths)
z_flux = zern_flux_zs(zeniths, azimuths)

cmap = mpl.colormaps['viridis']
bz1 = [np.min(z_fission), np.max(z_fission)]
bz2 = [np.min(z_flux), np.max(z_flux)]

nz1 = mpl.colors.Normalize(bz1[0], bz1[1] * 0.8)
nz2 = mpl.colors.Normalize(bz2[0], bz2[1])

smz1 = ScalarMappable(norm=nz1, cmap=cmap)
smz2 = ScalarMappable(norm=nz2, cmap=cmap)

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, r, z_fission, levels=1000)
plt.colorbar(smz1, ax=ax)
plt.savefig("zernike-kappa-fission.png", dpi=600)
plt.close()

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, r, z_flux, levels=1000)
plt.colorbar(smz2, ax=ax)
plt.savefig("zernike-flux.png", dpi=600)
plt.close()

# plotting polygon
poly_fission_df = poly_fission.get_pandas_dataframe()
poly_flux_df = poly_flux.get_pandas_dataframe()

poly_fission_n = poly_fission_df['mean']
poly_flux_n = poly_flux_df['mean']

poly_fission_zs = openmc.Zernike(poly_fission_n, radius = polygon_radius)
poly_flux_zs = openmc.Zernike(poly_flux_n, radius = polygon_radius)

num_sides = 6
alpha = np.pi / num_sides

def r_alpha(theta):
    drop = (theta + alpha) / (2 * alpha)
    u_alpha = theta - drop.astype(int) * (2 * alpha)
    return polygon_radius * np.cos(alpha) / np.cos(u_alpha)

var_radius = r_alpha(theta)
rp = r * var_radius / polygon_radius

k_fission = poly_fission_zs(zeniths, azimuths)
k_flux = poly_flux_zs(zeniths, azimuths)

bp1 = [np.min(k_fission), np.max(k_fission)]
bp2 = [np.min(k_flux), np.max(k_flux)]

np1 = mpl.colors.Normalize(bz1[0], bz1[1])
np2 = mpl.colors.Normalize(bz2[0], bz2[1])

smp1 = ScalarMappable(norm=nz1, cmap=cmap)
smp2 = ScalarMappable(norm=nz2, cmap=cmap)

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, rp, k_fission, levels=1000)
plt.colorbar(smp1, ax=ax)
plt.savefig("polygon-kappa-fission.png", dpi=600)
plt.close()

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, rp, k_flux, levels=1000)
plt.colorbar(smp2, ax=ax)
plt.savefig("polygon-flux.png", dpi=600)
plt.close()
