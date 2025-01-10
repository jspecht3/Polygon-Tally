import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib import patches as mp
from hex_rings import cs, pitch
from calculating_weights import values, max, min

# dimesions
radius = pitch / 2 / np.cos(np.radians(30))

# plotting
# bounds = [0.00298973, 0.0172874]
bounds = [min, max]
diff = bounds[1] - bounds[0]

cmap = mpl.colormaps['viridis']
norm = mpl.colors.Normalize(bounds[0], bounds[1])
sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])


def random():
    return diff * np.random.rand() + bounds[0]


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

plt.savefig("poly-flux-hex.png", dpi=600)
plt.show()
