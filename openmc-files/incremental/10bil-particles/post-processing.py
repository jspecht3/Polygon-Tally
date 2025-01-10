import openmc
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv

# importing
with open('transfer.csv', newline='') as file:
    reader = csv.reader(file, delimiter=' ')

    for row in reader:
        if row[0] == "mesh_dimensions":
            mesh_dimensions = (int(row[1]), int(row[2]))
        if row[0] == "pin_f2f":
            pin_cell_f2f = float(row[1])
        if row[0] == "polygon_radius":
            polygon_radius = float(row[1])

# statepoint
sp = openmc.StatePoint('statepoint.100000.h5')

mesh_tally = sp.tallies[1]
zern_tally = sp.tallies[2]
poly_tally = sp.tallies[3]

mesh_fission = mesh_tally.get_slice(scores=['kappa-fission'])
mesh_flux = mesh_tally.get_slice(scores=['flux'])

zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_flux = zern_tally.get_slice(scores=['flux'])

poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
poly_flux = poly_tally.get_slice(scores=['flux'])

mesh_fission.mean.shape = mesh_dimensions
mesh_flux.mean.shape = mesh_dimensions

# plotting
## mesh
plot = plt.subplot()
plot.imshow(np.flipud(mesh_fission.mean))
plt.savefig("mesh-kappa-fission.png", dpi=600)
plt.close()

plot = plt.subplot()
plot.imshow(np.flipud(mesh_flux.mean))
plt.savefig("mesh-flux.png", dpi=600)
plt.close()

## zernike
zern_fission_df = zern_fission.get_pandas_dataframe()
zern_flux_df = zern_flux.get_pandas_dataframe()

zern_fission_n = zern_fission_df['mean']
zern_flux_n = zern_flux_df['mean']

zern_fission_zs = openmc.Zernike(zern_fission_n, radius = polygon_radius)
zern_flux_zs = openmc.Zernike(zern_flux_n, radius = polygon_radius)

azimuths = np.linspace(0, 2*np.pi, 100)
zeniths = np.linspace(0, polygon_radius, 100)
r,theta = np.meshgrid(zeniths, azimuths)

z_fission = zern_fission_zs(zeniths, azimuths)
for i in range(len(z_fission)):
    for j in range(len(z_fission[i])):
        z_fission[i][j] /= 1e7
z_flux = zern_flux_zs(zeniths, azimuths)

### plotting
fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, r, z_fission, cmap="coolwarm", vmin=-0.005, vmax=0.021)
fig.colorbar(plot, pad=0.075)
plt.savefig("zernike-kappa-fission.png", dpi=600)
plt.close()

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, r, z_flux, cmap="coolwarm", vmin=-0.005, vmax=0.021)
fig.colorbar(plot, pad=0.075)
plt.savefig("zernike-flux.png", dpi=600)
plt.close()

## polygon
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
for i in range(len(k_fission)):
    for j in range(len(k_fission[i])):
        k_fission[i][j] /= 1e7
k_flux = poly_flux_zs(zeniths, azimuths)

### plotting
fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, rp, k_fission, cmap='coolwarm', vmin=-0.005, vmax=0.021)
fig.colorbar(plot, pad=0.075)
plt.savefig("polygon-kappa-fission.png", dpi=600)
plt.close()

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, rp, k_flux, cmap='coolwarm', vmin=-0.005, vmax=0.021)
fig.colorbar(plot, pad=0.075)
plt.savefig("polygon-flux.png", dpi=600)
plt.close()
