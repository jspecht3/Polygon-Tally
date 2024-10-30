import numpy as np
from polygon_tally.zernike_like import *
import matplotlib.pyplot as plt

# defining k basis
n, m = 0, 0
num_sides = 4
polygon_radius = 1
alpha = np.pi / num_sides
k = KBasis(n,m,num_sides,polygon_radius)

# coords
theta = np.linspace(0,2*np.pi,1000)
bleh = np.degrees(theta)

u_alpha = k.u_alpha_theta(theta)
r_alpha = k.variable_radius_theta(theta)

corners = []
labels = []
for i in range(4):
    corner = np.degrees(alpha + 2*i*alpha)
    corners.append(corner)
    labels.append(f"{corner}$\degree$")

plt.plot(bleh, u_alpha, label="U", ls=(0,(5,3)))
plt.plot(bleh, r_alpha, label="R")

plt.xlabel("Theta [${\degree}$]")
plt.xticks(ticks = corners, labels = labels)

plt.legend(framealpha = 1)
plt.savefig("ur-alpha.png", dpi=600)
plt.close()

# plotting the square
k._gen_all_coords(2500)
r,theta = k._r, k._theta

fig, ax = plt.subplots(subplot_kw=dict(projection="polar"))
plot = ax.contourf(theta,r,r*0)
plt.axvline(-alpha, color='k', ls=(0,(5,3)), label="First Sector")
plt.axvline(alpha, color='k', ls=(0,(5,3)))
plt.legend(loc=(0.75,0.95))
plt.savefig("square.png", dpi=600)
plt.close()
