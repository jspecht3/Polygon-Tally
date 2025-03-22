import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib import patches as mp
from hex_rings import cs, pitch
from max_error import (zern_flux_xy, zern_fission_xy, poly_flux_xy,
    poly_fission_xy)
from calc_max_err import calc_weights
from parse_csv import flux_vals, fission_vals
from functions import (zern_flux_err, zern_fission_err, poly_flux_err,
    poly_fission_err)

# dimensions
radius = pitch / 2 / np.cos(np.radians(30))


# plotting
def max_error(func, dx):
    tally_type, qoi, _ = func.__name__.split('_')

    if qoi == "flux":
        parsed = flux_vals.copy()
    if qoi == "fission":
        parsed = fission_vals.copy()
    elif qoi != "flux" and qoi != "fission":
        raise ValueError(f"{qoi}")

    # normalizing parsed data
    total_unc = 0
    for ring in cs:
        for cell_num in range(len(cs[ring])):
            total_unc += parsed[ring][cell_num]

    for ring in cs:
        for cell_num in range(len(cs[ring])):
            parsed[ring][cell_num] /= total_unc

    maxes = []
    # looping, momento moment
    for order in range(16):
        values, _, _ = calc_weights(func, dx, order)

        total_neph = 0
        for ring in cs:
            for cell_num in range(len(cs[ring])):
                total_neph += values[ring][cell_num]

        for ring in cs:
            for cell_num in range(len(cs[ring])):
                values[ring][cell_num] /= total_neph

        max_err = -1
        for ring in cs:
            for cell_num in range(len(cs[ring])):
                vle = values[ring][cell_num]
                psd = parsed[ring][cell_num]
                diff = abs(vle - psd)
                if diff > max_err:
                    max_err = diff

        maxes.append(max_err)

    print(f"{tally_type}-{qoi} : {maxes}")
    if tally_type == "zern":
        glip = "Zernike"
    if tally_type == "poly":
        glip = "Polygon"

    if qoi == "flux":
        glub = "Flux"
    if qoi == "fission":
        glub = "Kappa-Fission"

    if tally_type == "zern" and qoi == "fission":
        ls = "solid"
    if tally_type == "zern" and qoi == "flux":
        ls = "dotted"
    if tally_type == "poly" and qoi == "fission":
        ls = "dashed"
    if tally_type == "poly" and qoi == "flux":
        ls = "dashdot"

    ax1.semilogy(maxes, label=f"{glub}, {glip}")#, ls=ls, color='r')

fig, ax1 = plt.subplots()
#ax2 = plt.twinx(ax1)

def plot_all_maxes(dx):

    max_error(zern_fission_xy, dx)
    max_error(poly_fission_xy, dx)
    max_error(zern_flux_xy, dx)
    max_error(poly_flux_xy, dx)
    """
    ax2.semilogy(zern_flux_err, label = "RE: Zern Fl", marker='o', color='b')
    ax2.semilogy(zern_fission_err, label = "RE: Zern KF", marker='s', color='b')
    ax2.semilogy(poly_flux_err, label = "RE: Poly, FL", marker='^', color='b')
    ax2.semilogy(poly_fission_err, label = "RE: Poly, KF", marker='x', color='b')
    """
    ax1.legend(loc='upper right')
    #ax2.legend(loc='upper left')
    plt.ylim(5e-5, 1e-1)
    plt.ylabel("Maximum Error")
    plt.xlabel("FET Radial Order")

    plt.grid('both')
    plt.savefig("./plots/max-err.png", dpi=600)
    plt.show()


plot_all_maxes(0.3)
