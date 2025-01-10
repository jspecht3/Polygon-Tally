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

# zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
# zern_flux = zern_tally.get_slice(scores=['flux'])

# poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
# poly_flux = poly_tally.get_slice(scores=['flux'])


# zernike fission info
zern_fission = zern_tally.get_slice(scores=['kappa-fission'])
zern_fission_df = zern_fission.get_pandas_dataframe()
zern_fission_n = zern_fission_df['mean']
zern_fission_zs = openmc.Zernike(zern_fission_n, radius=1)

# zernike flux info
zern_flux = zern_tally.get_slice(scores=['flux'])
zern_flux_df = zern_flux.get_pandas_dataframe()
zern_flux_n = zern_flux_df['mean']
zern_flux_zs = openmc.Zernike(zern_flux_n, radius=1)

# polygon fission info
poly_fission = poly_tally.get_slice(scores=['kappa-fission'])
poly_fission_df = poly_fission.get_pandas_dataframe()
poly_fission_n = poly_fission_df['mean']
poly_fission_zs = openmc.Zernike(poly_fission_n, radius=1)

# polygon flux info
poly_flux = poly_tally.get_slice(scores=['flux'])
poly_flux_df = poly_flux.get_pandas_dataframe()
poly_flux_n = poly_flux_df['mean']
poly_flux_zs = openmc.Zernike(poly_flux_n, radius=1)


# functions for zern tally
def zern_fission_xy(x, y):
    r = (x**2 + y**2)**(1/2) / polygon_radius
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))
    return zern_fission_zs(r, theta) * np.pi


def zern_flux_xy(x, y):
    r = (x**2 + y**2)**(1/2) / polygon_radius
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))
    return zern_flux_zs(r, theta) * np.pi


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
    return poly_fission_zs(rho, theta) * np.pi


def poly_flux_xy(x, y):
    r = (x**2 + y**2)**(1/2)
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))

    var_radius = r_alpha(theta)
    rho = r / var_radius
    return poly_flux_zs(rho, theta) * np.pi


def integrate_circle(func, ms):
    r = np.linspace(0, polygon_radius, ms)
    theta = np.linspace(0, 2 * np.pi, ms)

    r, theta = np.meshgrid(r, theta)

    dr = r[:, 1]
    _, dr = np.meshgrid(dr, dr)
    dtheta = theta[1, 1]

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    integral = 0
    for xi, yi, ri, dri in zip(x, y, r, dr):
        integral += func(xi, yi) * ri * dri * dtheta
    print(f"Integral of circle {func.__name__} is {integral}")
    return


def integrate_hex(func, radius, ms):
    pitch = 2 * radius * np.cos(np.degrees(30))
    h = pitch / 2

    ys = np.linspace(-h, h, ms)
    n = len(ys)
    dy = y[1] - y[0]

    def get_xs(y):
        xr = radius - (radius * abs(y) / 2 / h)
        xl = - radius + (radius * abs(y) / 2 / h)

        m = abs(int(n * (1 - abs(y) / 2 / h)))
        xs = np.linspace(xl, xr, m)
        dx = xs[1] - xs[0]
        return xs, dx

        integral = 0
        for y in ys:
            xs, dx = gen_xs(y)
            for x in xs:
                integral += func(x, y) * dx * dy


def unity(x, y):
    return 1


# integrate_circle(unity, 10000)

integrate_circle(zern_flux_xy, 5000)
integrate_circle(zern_fission_xy, 5000)
