import openmc
import openmc.lib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
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
sp = openmc.StatePoint('statepoint.500.h5')

mesh_tally = sp.tallies[1]
zern_tally = sp.tallies[2]
poly_tally = sp.tallies[3]
dist_tally = sp.tallies[4]

mesh_fission = mesh_tally.get_slice(scores=['kappa-fission'])
mesh_flux = mesh_tally.get_slice(scores=['flux'])

zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_flux = zern_tally.get_slice(scores=['flux'])

poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
poly_flux = poly_tally.get_slice(scores=['flux'])

dist_fission = dist_tally.get_slice(scores=['kappa-fission'])
dist_flux = dist_tally.get_slice(scores=['flux'])


# plotting
## mesh
mesh_fission.mean.shape = mesh_dimensions
mesh_flux.mean.shape = mesh_dimensions

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
z_flux = zern_flux_zs(zeniths, azimuths)

### plotting
fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, r, z_fission)
plt.savefig("zernike-kappa-fission.png", dpi=600)
plt.close()

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, r, z_flux)
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
k_flux = poly_flux_zs(zeniths, azimuths)

### plotting
fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, rp, k_fission)
plt.savefig("polygon-kappa-fission.png", dpi=600)
plt.close()

fig, ax = plt.subplots(subplot_kw = dict(projection="polar"))
plot = ax.contourf(theta, rp, k_flux)
plt.savefig("polygon-flux.png", dpi=600)
plt.close()

'''
# distribcell
pellet_diameter = 5.5685 / 10
pin_pitch = 1.180 * pellet_diameter
pin_radius = pin_pitch / 3**(0.5)

dist_df = dist_fission.get_pandas_dataframe()

x = dist_df['level 3']['lat']['x']
y = dist_df['level 3']['lat']['y']
means = dist_df['mean']

fig, ax = plt.subplots()
ax.set_aspect('equal')

for i,j in zip(x,y):
    if j % 2 == 0:
        i += 1/2
#        j += 1 / 3**(1/2)
    hex = RegularPolygon(
            (i,j),
            numVertices = 6,
            radius = 3**(-1/2),
            alpha = 0.2,
            edgecolor = 'k'
    )
    ax.add_patch(hex)
plt.autoscale(enable = True)
plt.savefig("hex.png", dpi=600)
plt.close()


resolution = (600, 600)
img = np.full(resolution, np.nan)
xmin, xmax = -12, 12
ymin, ymax = -12, 12
with openmc.lib.run_in_memory():

    for row, y in enumerate(np.linspace(ymin, ymax, resolution[0])):
        for col, y in enumerate(np.linspace(xmin, xmax, resolution[1])):
            try:
                cell, distrib_cell_index = openmc.lib.find_cell((x,y,0))
            except openmc.exceptions.GeometryError:
                continute

            print(cell.id)
            if cell.id == 2:
                img[row, col] = fission[distribcell_index]

options = {
        'origin': 'lower',
        'extent': (xmin, xmax, ymin, ymax)
}

plt.imshow(img, **options)
plt.colorbar()
plt.savefig("distrib-fission.png", dpi=600)
'''
