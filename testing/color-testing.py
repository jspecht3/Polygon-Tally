import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mp
from matplotlib.cm import ScalarMappable
import numpy as np

cmap = mpl.colormaps['viridis']

xs = np.linspace(0, 1, 15)
ys = np.linspace(0, 0, 15)
data = np.linspace(0.23, .74, 15)

norm = mpl.colors.Normalize(min(data), max(data))

sm = ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

fig, ax = plt.subplots()

for i in range(len(xs)):
    color = cmap(norm(data[i]))
    ax.scatter(xs[i], ys[i], color=color)

plt.colorbar(sm, ax=ax)
plt.show()
