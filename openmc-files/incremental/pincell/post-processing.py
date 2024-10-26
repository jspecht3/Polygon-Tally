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

# statepoint
sp = openmc.StatePoint('statepoint.100.h5')

mesh_tally = sp.tallies[1]
zern_tally = sp.tallies[2]

mesh_fission = mesh_tally.get_slice(scores=['kappa-fission'])
mesh_flux = mesh_tally.get_slice(scores=['flux'])

zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_flux = zern_tally.get_slice(scores=['flux'])

mesh_fission.mean.shape = mesh_dimensions
mesh_flux.mean.shape = mesh_dimensions

# plotting
## mesh
plot = plt.subplot()
plot.imshow(mesh_fission.mean)
plt.savefig("mesh-kappa-fission.png", dpi=600)
plt.close()

plot = plt.subplot()
plot.imshow(mesh_flux.mean)
plt.savefig("mesh-flux.png", dpi=600)
plt.close()

## zernike
zern_fission_df = zern_fission.get_pandas_dataframe()
zern_flux_df = zern_flux.get_pandas_dataframe()

zern_fission_n = zern_fission_df['mean']
zern_flux_n = zern_flux_df['mean']

zern_fission_zs = openmc.Zernike(
        zern_fission_n,
        radius = pin_cell_f2f / 2)
zern_flux_zs = openmc.Zernike(
        zern_flux_n,
        radius = pin_cell_f2f / 2)

azimuths = np.linspace(0, 2*np.pi, 100)
zeniths = np.linspace(0, pin_cell_f2f / 2, 100)
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
