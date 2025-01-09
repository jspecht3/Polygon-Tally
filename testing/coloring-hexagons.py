import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon as rp
import numpy as np

polygon = rp((0, 0), 6, radius=3, orientation=np.pi/2)

fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.add_artist(polygon)
ax.set_axis_off()

plt.xlim(-5, 5)
plt.ylim(-5, 5)

plt.show()
