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


# getting zernike fission function
zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_fission_df = zern_fission.get_pandas_dataframe()
zern_fission_n = zern_fission_df['mean']
zern_fission_zs = openmc.Zernike(zern_fission_n, radius=polygon_radius)

print(zern_fission_n[0])
