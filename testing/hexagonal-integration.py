import numpy as np
import matplotlib.pyplot as plt

# want a function that takes a fn(x,y), a center, and a radius
# numerically integrates the value over the hexagon w/ uniform mesh spacing

def get_ys(n, h, b): return np.linspace(b - h, b + h, n)


def get_xs(r, y, n, h, a):
    xr = r + a - (r / 2 / h) * abs(y)
    xl = 2 * a - xr

    m = int(n * (1 - abs(y) / 2 / h))
    return np.linspace(xl, xr, m)


def get_coords(r, dx, a = 0, b = 0):
    """generates a roughly evenly spaced quadrature set over a hexagon

    Parameters
    ----------
    r : float
        radius of the hexagon
    dx : float
        quadrature set spacing
    a : float
        x coordinate of the hexagon center
    b : float
        y coordiante of the hexagon center
    """

    h = r / 2 * np.tan(np.radians(60)) # height of hexagon
    n = int(2 * h / dx) # max number of quad points, quad points at y = b

    ys = get_ys(n, h, b)
    points = []

    for y in ys:
        xs = get_xs(r, y, n, h, a)
        for x in xs:
            point = (x, y)
            points.append(point)

    return points


coords = get_coords(5, 1)
for coord in coords:
    x,y = coord
    plt.scatter(x,y)
plt.show()
