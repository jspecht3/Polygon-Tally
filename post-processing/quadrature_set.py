import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon as rp
from hex_rings import pitch, radius, cs

# dimesions
h = pitch / 2


def get_ys(n, y0):
    return np.linspace(y0 - h, y0 + h, n)


def get_xs(n, x0, y0, y):
    xr = x0 + radius - (radius * abs(y - y0) / 2 / h)
    xl = x0 - radius + (radius * abs(y - y0) / 2 / h)

    m = abs(int(n * (1 - abs(y - y0) / 2 / h)))
    return np.linspace(xl, xr, m)


def quad_set(center, dx):
    x0, y0 = center
    n = int(2 * h / dx)

    ys = get_ys(n, y0)
    points = []

    count = 0

    for y in ys:
        xs = get_xs(n, x0, y0, y)
        count += len(xs)
        for x in xs:
            points.append(np.array([x, y]))

    return points, count


# plotting
def plotter():
    center = np.array([0, 10])

    fig, ax = plt.subplots()
    ax.set_aspect('equal')
    ax.add_artist(rp((center[0], center[1]), 6, radius=radius, orientation=np.pi/2, 
                     color='b', ec='r'))

    coords, count = quad_set(center, 0.1)
    print(count)

    for coord in coords:
        x, y = coord
        plt.scatter(x, y, color='k')

    plt.savefig("quadrature-set.png", dpi=600)
    plt.show()

# plotter()
# print(radius, h)