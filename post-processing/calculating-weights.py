from functions import zern_fission_xy, poly_fission_xy, polygon_radius
from quadrature_set import quad_set, plotter
from hex_rings import cs, radius

values = {}
num_quad = {}

for ring in cs:
    for cell_num in range(len(cs[ring])):
        cell_id = f"cell_{ring}_{cell_num}"
        center = cs[ring][cell_num]

        quad_points, count = quad_set(center, 0.1)
        value = 0

        for quad_point in quad_points:
            print(quad_point)
