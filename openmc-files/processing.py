import openmc
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# importing the filters
from openmc.filter import MeshFilter, DistribcellFilter
from openmc.filter_expansion import ZernikeFilter

# post-processing
## getting scores
sp = openmc.StatePoint('statepoint.10000.h5')

for tally_key in sp.tallies:
    tally = sp.tallies[tally_key]
    filter_type = type(tally.filters[0])
    filter_ = tally.filters[0]

    if filter_type == MeshFilter and tally.scores == ['kappa-fission']:
        mesh_heat_scores = tally
    if filter_type == MeshFilter and tally.scores == ['flux']:
        mesh_flux_scores = tally 

    if filter_type == DistribcellFilter and tally.scores == ['kappa-fission']:
        distrib_heat_scores = tally
    if filter_type == DistribcellFilter and tally.scores == ['flux']:
        distrib_flux_scores = tally

    if filter_type == ZernikeFilter and tally.scores == ['kappa-fission'] and filter_.num_sides == 0:
        zernike_heat_scores = tally
    if filter_type == ZernikeFilter and tally.scores == ['flux'] and filter_.num_sides == 0:
        zernike_flux_scores = tally

    if filter_type == ZernikeFilter and tally.scores == ['kappa-fission'] and filter_.num_sides == 6:
        polygon_heat_scores = tally
    if filter_type == ZernikeFilter and tally.scores == ['flux'] and filter_.num_sides == 6:
        polygon_flux_scores = tally

## mesh
mesh_len = 40
mesh_dimensions = (mesh_len, mesh_len)

mesh_heat_scores.std_dev.shape = mesh_dimensions
mesh_heat_scores.mean.shape = mesh_dimensions

mesh_flux_scores.std_dev.shape = mesh_dimensions
mesh_flux_scores.mean.shape = mesh_dimensions

mesh_heat_plot = plt.subplot()
mesh_heat_plot.imshow(mesh_heat_scores.mean)
plt.savefig("post/mesh_heat.png", dpi=600)
plt.close()

mesh_flux_plot = plt.subplot()
mesh_flux_plot.imshow(mesh_flux_scores.mean)
plt.savefig("post/mesh_flux.png", dpi=600)
plt.close()


## distribcell is something to as novak about


## zernike params
pellet_diameter = 5.5685 / 10
pin_pitch = 1.180 * pellet_diameter
polygon_radius = pin_pitch * 3**(0.5) * 10

## circular zernike
zernike_heat_df = zernike_heat_scores.get_pandas_dataframe()
zernike_flux_df = zernike_flux_scores.get_pandas_dataframe()

zernike_heat_n = zernike_flux_df['mean']
zernike_flux_n = zernike_flux_df['mean']

zernike_heat_zernikes = openmc.Zernike(zernike_heat_n, radius = polygon_radius)
zernike_flux_zernikes = openmc.Zernike(zernike_flux_n, radius = polygon_radius)

### coordinates
azimuths = np.radians(np.linspace(0,360,100))
zeniths = np.linspace(0,polygon_radius,100)
r,theta = np.meshgrid(zeniths, azimuths)

z_heat_values = zernike_heat_zernikes(zeniths, azimuths)
z_flux_values = zernike_flux_zernikes(zeniths, azimuths)

### plotting
fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
plot = ax.contourf(theta,r,z_heat_values,cmap="jet")
cbar = fig.colorbar(plot)
plt.savefig("post/z_heat.png", dpi=600)
plt.close()

fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
plot = ax.contourf(theta,r,z_flux_values,cmap="jet")
cbar = fig.colorbar(plot)
plt.savefig("post/z_flux.png", dpi=600)
plt.close()


## polygon zernike
### parameters
num_sides = 6
alpha = np.pi / num_sides

### getting data
polygon_heat_df = polygon_heat_scores.get_pandas_dataframe()
polygon_flux_df = polygon_flux_scores.get_pandas_dataframe()

polygon_heat_n = polygon_flux_df['mean']
polygon_flux_n = polygon_flux_df['mean']

polygon_heat_zernikes = openmc.Zernike(polygon_heat_n, radius = polygon_radius)
polygon_flux_zernikes = openmc.Zernike(polygon_flux_n, radius = polygon_radius)

### transforming
def r_alpha(theta):
    drop = (theta + alpha) / (2 * alpha)
    u_alpha = theta - drop.astype(int) * (2 * alpha) 

    return polygon_radius * np.cos(alpha) / np.cos(u_alpha)

var_radius = r_alpha(theta)
rp = r * var_radius / polygon_radius

### values
k_heat_values = polygon_heat_zernikes(zeniths, azimuths)
k_flux_values = polygon_flux_zernikes(zeniths, azimuths)

### plotting
fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
plot = ax.contourf(theta,rp,k_heat_values,cmap="jet")
cbar = fig.colorbar(plot)
plt.savefig("post/k_heat.png")
plt.close()

fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
plot = ax.contourf(theta,rp,k_flux_values,cmap="jet")
cbar = fig.colorbar(plot)
plt.savefig("post/k_flux.png")
plt.close()
