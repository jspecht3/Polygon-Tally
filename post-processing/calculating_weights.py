from functions import (zern_fission_xy, zern_flux_xy, poly_fission_xy,
poly_flux_xy, polygon_radius)
from quadrature_set import quad_set, plotter
from hex_rings import cs, radius
import numpy as np
import time


def calc_weights(func, dx):
    tally_type, qoi, _ = func.__name__.split('_')

    if qoi == "fission":
        normalizer = 12826935.117376836
    if qoi == "flux":
        normalizer = 2.0552252049066863e-12
    if qoi != "fission" and qoi != "flux":
        raise ValueError("put the fission or flux in the bag lil bro")

    t0 = time.time()
    values = {}
    num_quad = 0

    area = 6 * radius**2 * np.cos(np.pi / 6) * np.sin(np.pi / 6)

    for ring in cs:
        values[ring] = []
        for cell_num in range(len(cs[ring])):
            cell_id = f"cell_{ring}_{cell_num}"
            center = cs[ring][cell_num]

            quad_points, count = quad_set(center, dx)
            value = 0

            for quad_point in quad_points:
                value += func(quad_point[0], quad_point[1])

            values[ring].append(value / count / area / normalizer)
            num_quad = count

    max = 1e-100
    min = 1e100
    # loc = np.array([0, 0])
    # mloc = np.array([0, 0])
    for ring in values:
        for cell in range(len(values[ring])):
            if values[ring][cell] > max:
                max = values[ring][cell]
                # loc = np.array([ring, cell])
            if values[ring][cell] < min:
                min = values[ring][cell]
                # mloc = np.array([ring, cell])

    # print(f"Max: {max}, Min: {min}")
    # print(min, mloc)
    # print([9, 18], values[9][18])
    # print(num_quad)

    # print(values)

    print(f"Quad Points: {num_quad}")
    tf = time.time()
    print(f"Finished calculating the {tally_type} {qoi} values in {tf - t0}s")

    return values, min, max


def area_calc(radius):
    return 6 * radius**2 * np.cos(np.pi / 6) * np.sin(np.pi / 6)


def print_areas():
    print(f"Little Hex Area: {area_calc(radius)}")
    print(f"Big Hex Area: {area_calc(polygon_radius)}")
    print(f"Big Circle Area: {np.pi * polygon_radius**2}")


# print_areas()
