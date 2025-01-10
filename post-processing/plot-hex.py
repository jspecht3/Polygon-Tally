import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import ScalarMappable
from matplotlib import patches as mp
from hex_rings import cs, pitch

# dimesions
radius = pitch / 2 / np.cos(np.radians(30))

# plotting
bounds = [0.00298973, 0.0172874]
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
    for cell in cs[ring]:
        x, y = cell
        color = cmap(norm(random()))
        polygon = mp.RegularPolygon((x, y), 6, radius=radius,
                                    orientation=np.pi/2, color=color, ec='r')
        ax.add_artist(polygon)

plt.xlim(-11, 11)
plt.ylim(-12, 12)
plt.colorbar(sm, ax=ax)

plt.savefig("plot-hex.png", dpi=600)
plt.show()
