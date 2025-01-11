import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib import patches as mp
from hex_rings import cs, pitch
from functions import (zern_flux_xy, zern_fission_xy, poly_flux_xy,
    poly_fission_xy)
from calculating_weights import calc_weights
from parse_csv import flux_vals, fission_vals

# dimensions
radius = pitch / 2 / np.cos(np.radians(30))


# plotting
def plot_func(func, dx):
    values, min, max = calc_weights(func, dx)
    tally_type, qoi, _ = func.__name__.split('_')

    total = 0
    for ring in cs:
        for cell_num in range(len(cs[ring])):
            x, y = cs[ring][cell_num]
            # color = cmap(norm(random()))
            total += values[ring][cell_num]

    # bounds = [min / total, max / total]
    # diff = bounds[1] - bounds[0]

    if qoi == "fission":
        bounds = [0, 0.018]
    if qoi == "flux":
        bounds = [0, 0.009]

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    cmap = mpl.colormaps['viridis']
    norm = mpl.colors.Normalize(bounds[0], bounds[1])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    for ring in cs:
        for cell_num in range(len(cs[ring])):
            x, y = cs[ring][cell_num]
            # color = cmap(norm(random()))
            weight = values[ring][cell_num] / total
            color = cmap(norm(weight))
            polygon = mp.RegularPolygon((x, y), 6, radius=radius,
                                        orientation=np.pi/2, color=color)
            ax.add_artist(polygon)

    plt.xlim(-11, 11)
    plt.ylim(-12, 12)
    plt.colorbar(sm, ax=ax)

    plt.savefig(f"plots/{tally_type}-{qoi}-hex.png", dpi=600)
    plt.close()


# plotting csv
def plot_csv(qoi):
    if qoi == "flux":
        parsed = flux_vals.copy()
    if qoi == "fission":
        parsed = fission_vals.copy()

    max = -1e10
    min = 1e10
    total = 0
    for ring in parsed:
        for cell_num in range(len(parsed[ring])):
            total += parsed[ring][cell_num]

    for ring in parsed:
        for cell_num in range(len(parsed[ring])):
            weight = parsed[ring][cell_num]
            if weight / total > max:
                max = weight / total
            if weight / total < min:
                min = weight / total

#    if qoi == 'fission':
#        bounds = [0, 0.005]
#    if qoi == 'flux':
#        bounds = [0, 3e9]
#    bounds = [min, max]

    if qoi == "fission":
        bounds = [0, 0.018]
    if qoi == "flux":
        bounds = [0, 0.009]

    cmap = mpl.colormaps['viridis']
    norm = mpl.colors.Normalize(bounds[0], bounds[1])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    for ring in cs:
        for cell_num in range(len(cs[ring])):
            x, y = cs[ring][cell_num]
            weight = parsed[ring][cell_num] / total
            color = cmap(norm(weight))
            polygon = mp.RegularPolygon((x, y), 6, radius=radius,
                                        orientation=np.pi/2, color=color)
            ax.add_artist(polygon)

    plt.xlim(-11, 11)
    plt.ylim(-12, 12)
    plt.colorbar(sm, ax=ax)

    plt.savefig(f"plots/{qoi}-csv.png", dpi=600)
    plt.close()


def plot_diffs(func, dx):
    values, min, max = calc_weights(func, dx)
    tally_type, qoi, _ = func.__name__.split('_')

    if qoi == "flux":
        parsed = flux_vals.copy()
        bounds = [0, 0.001]
    if qoi == "fission":
        parsed = fission_vals.copy()
        bounds = [0, 0.005]
    elif qoi != "flux" and qoi != "fission":
        raise ValueError(f"{qoi}")

    total1 = 0
    total2 = 0
    for ring in cs:
        for cell_num in range(len(cs[ring])):
            total1 += values[ring][cell_num]
            total2 += parsed[ring][cell_num]

    for ring in cs:
        for cell_num in range(len(cs[ring])):
            values[ring][cell_num] /= total1
            parsed[ring][cell_num] /= total2

    diff = {}
    for ring in parsed:
        diff[ring] = []
        for cell_num in range(len(parsed[ring])):
            vle = values[ring][cell_num]
            psd = parsed[ring][cell_num]
            diff[ring].append(abs(vle - psd))

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    cmap = mpl.colormaps['viridis']
    norm = mpl.colors.Normalize(bounds[0], bounds[1])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    max_diff = 0
    avg_diff = 0
    for ring in cs:
        for cell_num in range(len(cs[ring])):
            x, y = cs[ring][cell_num]
            weight = diff[ring][cell_num]
            avg_diff += weight
            if weight > max_diff:
                max_diff = weight
            color = cmap(norm(weight))
            polygon = mp.RegularPolygon((x, y), 6, radius=radius,
                                        orientation=np.pi/2, color=color)
            ax.add_artist(polygon)

    plt.xlim(-11, 11)
    plt.ylim(-12, 12)
    plt.colorbar(sm, ax=ax)

    plt.savefig(f"plots/diffs/{tally_type}-{qoi}.png", dpi=600)
    plt.close()

    print(f"----- {qoi} {tally_type} -----")
    print(f"max difference : {max_diff}")
    print(f"avg difference : {avg_diff / 271}")


def plot_all_func(dx):
    plot_func(zern_fission_xy, dx)
    plot_func(zern_flux_xy, dx)
    plot_func(poly_fission_xy, dx)
    plot_func(poly_flux_xy, dx)


def plot_all_csv():
    plot_csv('fission')
    plot_csv('flux')


def plot_all_diffs(dx):
    plot_diffs(zern_fission_xy, dx)
    plot_diffs(zern_flux_xy, dx)
    plot_diffs(poly_fission_xy, dx)
    plot_diffs(poly_flux_xy, dx)


# plot_all_func(0.04)
# plot_all_csv()
plot_all_diffs(0.04)
