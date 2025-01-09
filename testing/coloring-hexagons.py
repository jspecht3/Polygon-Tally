import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon as rp
import numpy as np

# dimensions
pellet_diameter = 5.5685 / 10
pitch = 1.180 * pellet_diameter
r = pitch / 2 / np.cos(np.radians(30))

polygon = rp((0, 0), 6, radius=r, orientation=np.pi/2, color='r', ec='b')

fig, ax = plt.subplots()
ax.set_aspect('equal')

ax.add_artist(polygon)

center = np.matrix([0, pitch])
for i in range(6):
    polygon = rp((center[0, 0], center[0, 1]), 6,
                 radius=r, orientation=np.pi/2, color='r', ec='b')
    ax.add_artist(polygon)
#    plt.scatter(center[0][0], center[0][1])

    print(center[0, 0])
    center *= np.matrix([
        [np.cos(np.degrees(60)), -np.sin(np.degrees(60))],
        [np.sin(np.degrees(60)), np.cos(np.degrees(60))]])

ax.set_axis_off()

plt.xlim(-2, 2)
plt.ylim(-2, 2)

plt.show()
