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
    return np.array([zern_fission_zs(r, theta)], dtype=float) * np.pi


def zern_flux_xy(x, y):
    r = (x**2 + y**2)**(1/2) / polygon_radius
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))
    return np.array([zern_flux_zs(r, theta)], dtype=float) * np.pi


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
    return np.array([poly_fission_zs(rho, theta)], dtype=float) * np.pi


def poly_flux_xy(x, y):
    r = (x**2 + y**2)**(1/2)
    theta = np.arctan2(y, x) + (2 * np.pi * (y < 0))

    var_radius = r_alpha(theta)
    rho = r / var_radius
    return np.array([poly_flux_zs(rho, theta)], dtype=float) * np.pi


def integrate_circle(func, ms):
    r = np.linspace(0, polygon_radius, ms)
    theta = np.linspace(0, 2 * np.pi, ms)

    r, theta = np.meshgrid(r, theta)

    dr = r[:, 1]
    _, dr = np.meshgrid(dr, dr)
    dtheta = theta[1, 1]

    x = r * np.cos(theta)
    y = r * np.sin(theta)

#    integral = np.sum(func(x, y) * r * dr * dtheta)
    integral = 0
    for xi, yi, ri, dri in zip(x, y, r, dr):
        integral += func(xi, yi) * ri * dri * dtheta
    print(f"Integral of circle {func.__name__} is {np.sum(integral)}")
    return


def integrate_hex(func, radius, ms):
    pitch = 2 * radius * np.cos(np.radians(30))
    h = pitch / 2
    q = radius * np.sin(np.radians(30))

    xs = np.linspace(-h, h, ms)
    dx = xs[1] - xs[0]

    slope = (radius - q) / h
    integral = 0

    for x in xs:
        if x <= 0:
            ymax = radius + slope * x
        if x > 0:
            ymax = radius - slope * x

        ys = np.linspace(-ymax, ymax, ms)
        dy = ys[1] - ys[0]

        #print(dx, dx, x, ys)

        for y in ys:
            integral += func(x, y) * dx * dy

    print(f"Hex {func.__name__}: {np.sum(integral)}")
    return


def unity(x, y):
    return 1


# integrate_circle(unity, 10000)

# integrate_circle(zern_flux_xy, 100)
# integrate_circle(zern_fission_xy, 100)

# integrate_hex(unity, radius, 500)
# integrate_hex(unity, polygon_radius, 500)

# integrate_hex(poly_flux_xy, polygon_radius, 100)
# integrate_hex(poly_fission_xy, polygon_radius, 100)


def fail():
    rho = np.linspace(0, 1, 100)
    phi = np.linspace(0, 1, 100)
    rho, phi = np.meshgrid(rho, phi)

    var_rad = r_alpha(phi)
    r = rho * var_rad
    theta = phi

    x = r * np.cos(theta)
    y = r * np.sin(theta)

    f = np.ones_like(x) * 100

    fig, ax = plt.subplots()
    # ax.contour(x, y, f)
    for xi, yi in zip(x, y):
        plt.scatter(x, y)

    plt.show()

    print(x, y, f)


# pitch = 2 * radius * np.cos(np.radians(30))
#h = pitch / 2
#q = radius * np.sin(np.radians(30))
#ms = 100

#xs = np.linspace(-h 
