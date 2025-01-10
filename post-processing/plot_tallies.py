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

    bounds = [min, max]
    # diff = bounds[1] - bounds[0]

    cmap = mpl.colormaps['viridis']
    norm = mpl.colors.Normalize(bounds[0], bounds[1])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    for ring in cs:
        for cell_num in range(len(cs[ring])):
            x, y = cs[ring][cell_num]
            # color = cmap(norm(random()))
            weight = values[ring][cell_num]
            color = cmap(norm(weight))
            polygon = mp.RegularPolygon((x, y), 6, radius=radius,
                                        orientation=np.pi/2, color=color)
            ax.add_artist(polygon)

    plt.xlim(-11, 11)
    plt.ylim(-12, 12)
    plt.colorbar(sm, ax=ax)

    plt.savefig(f"plots/{tally_type}-{qoi}-hex.png", dpi=600)
    plt.close()


# plotting the difference
def plot_diff_abs(func, dx):
    values, _, _ = calc_weights(func, dx)
    tally_type, qoi, _ = func.__name__.split('_')

    if qoi == "flux":
        parsed = flux_vals
    if qoi == "fission":
        parsed = fission_vals

    diffs = {}
    max = -1e10
    min = 1e10
    for ring in values:
        diffs[ring] = []
        for cell_num in range(len(values[ring])):
            diff = abs(values[ring][cell_num] - parsed[ring][cell_num]) 
            if diff > max:
                max = diff
            if diff < min:
                min = diff
            diffs[ring].append(diff)

    if qoi == 'fission':
        bounds = [0, 0.005]
    if qoi == 'flux':
        bounds = [0, 3e9]
#    bounds = [min, max]

    cmap = mpl.colormaps['viridis']
    norm = mpl.colors.Normalize(bounds[0], bounds[1])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    for ring in cs:
        for cell_num in range(len(cs[ring])):
            x, y = cs[ring][cell_num]
            weight = diffs[ring][cell_num]
            color = cmap(norm(weight))
            polygon = mp.RegularPolygon((x, y), 6, radius=radius,
                                        orientation=np.pi/2, color=color)
            ax.add_artist(polygon)

    plt.xlim(-11, 11)
    plt.ylim(-12, 12)
    plt.colorbar(sm, ax=ax)

    plt.savefig(f"plots/abs/{tally_type}-{qoi}-abs-diff.png", dpi=600)
    plt.close()


def plot_diff_rel(func, dx):
    values, _, _ = calc_weights(func, dx)
    tally_type, qoi, _ = func.__name__.split('_')

    if qoi == "flux":
        parsed = flux_vals
    if qoi == "fission":
        parsed = fission_vals

    diffs = {}
    max = -1e10
    min = 1e10
    for ring in values:
        diffs[ring] = []
        for cell_num in range(len(values[ring])):
            diff = abs((values[ring][cell_num] - parsed[ring][cell_num]) / parsed[ring][cell_num])
            if diff > max:
                max = diff
            if diff < min:
                min = diff
            diffs[ring].append(diff)

    if qoi == 'fission':
        bounds = [0, 1.4]
    if qoi == 'flux':
        bounds = [0, 0.25]
#    bounds = [min, max]

    cmap = mpl.colormaps['viridis']
    norm = mpl.colors.Normalize(bounds[0], bounds[1])
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    fig, ax = plt.subplots()
    ax.set_aspect('equal')

    for ring in cs:
        for cell_num in range(len(cs[ring])):
            x, y = cs[ring][cell_num]
            weight = diffs[ring][cell_num]
            color = cmap(norm(weight))
            polygon = mp.RegularPolygon((x, y), 6, radius=radius,
                                        orientation=np.pi/2, color=color)
            ax.add_artist(polygon)

    plt.xlim(-11, 11)
    plt.ylim(-12, 12)
    plt.colorbar(sm, ax=ax)

    plt.savefig(f"plots/rel/{tally_type}-{qoi}-rel-diff.png", dpi=600)
    plt.close()


def plot_all_func(dx):
    plot_func(zern_fission_xy, dx)
    plot_func(zern_flux_xy, dx)
    plot_func(poly_fission_xy, dx)
    plot_func(poly_flux_xy, dx)


def plot_all_diff_abs(dx):
    plot_diff_abs(zern_fission_xy, dx)
    plot_diff_abs(zern_flux_xy, dx)
    plot_diff_abs(poly_fission_xy, dx)
    plot_diff_abs(poly_flux_xy, dx)


def plot_all_diff_rel(dx):
#    plot_diff_rel(zern_fission_xy, dx)
    plot_diff_rel(zern_flux_xy, dx)
#    plot_diff_rel(poly_fission_xy, dx)
    plot_diff_rel(poly_flux_xy, dx)


# plot_all_func(0.03)
# plot_all_diff_abs(0.04)
plot_all_diff_rel(0.04)
